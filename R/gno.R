## Genetype Management

#' minor allele frequencies
#'
#' @param g genotype matrix, N row samples and M column variants
#' @return MAF of M variants
maf <- function(g) {a <- colMeans(g, na.rm=TRUE) / 2; pmin(a, 1 - a)}

#' allele sandard deviation
#'
#' @param g genotype matrix, N row samples and M column variants
#' @return SD of M variants
asd <- function(g) apply(g, 2, sd)

#' Cached Real Genotype
#' 
#' A segment taken from from 1000 genome project.
#'
#' A continuous segment from chromosome 17  region q12, with no missed calls and
#' minor allele frequency no less than 0.05.
#'
#' @name c17q12
#' @noRd
if(!exists("c17"))
{
    c17 <- readRDS(system.file("extdata", '17q12.rds', package="simgen"))
    n17 <- nrow(c17)
    m17 <- ncol(c17)
    f17 <- maf(c17)
}

#' Genotype from 1000 Genome Project
#'
#' Randomly draw a segment from the 1000 Genome Project.
#'
#' @param N number of samples
#' @param L number of SNPs
#' @param kgp.con draw variants continousely (def=1)
#' @param kgp.psd positive definite threshold
#' @param kgp.ucr upper correlation threshold
#' @param kgp.maf lower minor allele frequency threshold
kgp <- function(N, L=1, kgp.con=1, kgp.psd=NULL, kgp.ucr=NULL, min.maf=NULL, drop=TRUE, ...)
{
    psd <- kgp.psd %||% sqrt(.Machine$double.eps)
    ucr <- kgp.ucr %||% 0.99
    MAF <- min.maf %||% 0.05

    ## L variants on demand, reserve 2 * L
    P <- min(L * 5, m17)
    while(TRUE)
    {
        i <- sample.int(n17, N, N > n17)          # N
        if(kgp.con)                               # P
            j <- seq(sample(m17 - P, 1) + 1, l=P) #
        else                                      #
            j <- sort(sample(m17, P))             #
        gmx <- c17[i, j, drop=FALSE]              # N x P
        gmx <- gmx[, maf(gmx) >= MAF, drop=FALSE] # min MAF
        if(NCOL(gmx) < L)
        {
            P <- min(P + L - NCOL(gmx), m17)
            cat("P = ", P, "\n", sep="")
            next
        }
        ldm <- cor(gmx)                             # P x P
        
        ## drop  SNP  to  enforce  non-linearity;  if there  are  less  than  L
        ## remaining, try again with a bigger reserve.
        gmx <- gmx[, ncv(ldm, ucr=ucr, ...), drop=FALSE]
        if(NCOL(gmx) < L)
        {
            P <- min(P + L - NCOL(gmx), m17)
            cat("P = ", P, "\n", sep="")
            next
        }

        ## select L variants now
        j <- seq(sample(ncol(gmx) - L, 1) + 1, l=L)
        gmx <- gmx[, j, drop=drop]
        ldm <- if(L > 1) cor(gmx) else 1

        ## if min(eigenvalue) < (threshold) * max(eigenvalue), try again
        egv <- try(eigen(ldm, TRUE, TRUE)$values)
        if (egv[L] < psd * egv[1L])
        {
            cat("Non-PSD!\n")
            next
        }
        break
    }
    gmx
}

#' Binomial Genotype
#'
#' @param N sample size
#' @param L number of variants
#' @param MAF minor allele frequency
#' @param drop return a vector instead of a matrix when L==1.
#' @param std standardize genotype variance to 1
#' @param ... additional arguments
#' @return N*L genotype matrix generated from two trail binomial 
bng <- function(N, L=1, MAF=NULL, drop=TRUE, std=FALSE, ...)
{
    MAF <- MAF %||% 0.05
    gmx <- rbinom(N * L, 2L, MAF)
    if(std)
        gmx <- gmx / sqrt(2 * MAF * (1 - MAF))
    if(L > 1 || !drop)
        dim(gmx) <- c(N, L)
    gmx
}

#' Discretizing dosage by Hardy-Weinberg Equilibrium
#'
#' @param x matrix of allele dosages in continuous scale
#' @param m vector of allele frequencies
#' @return allele dosage in  {0, 1, 2}  with frequency {m^2, 2*m*(1-m), (1-m)^2}
#' @noRd
as.genotype <- function(x, m=NULL, ...)
{
    if (is.null(m))
    {
        m <- f17[sample(m17 - ncol(x), 1) + seq(ncol(x))]
    }
    for(j in seq(ncol(x)))
    {
        m2 <- m[j] * m[j]           # HWE freq of 0 allele
        mn <- 2 * m[j] * (1 - m[j]) # HWE freq of 1 allele
        qt <- quantile(x[, j], c(m2, m2 + mn))
        x[, j] <- 2 - (x[, j] > qt[1]) - (x[, j] > qt[2])
    }
    x
}
