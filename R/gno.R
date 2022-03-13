## Genetype Management

#' minor allele frequencies
#'
#' @param g genotype matrix, N row samples and M column variants
#' @return MAF of M variants
#' @export
maf <- function(g) {a <- colMeans(g, na.rm=TRUE) / 2; pmin(a, 1 - a)}

#' allele sandard deviation
#'
#' @param g genotype matrix, N row samples and M column variants
#' @return SD of M variants
#' @export
asd <- function(g) apply(g, 2, stats::sd)

#' retain non-colinear variables
#'
#' Trim a correlation  matrix so the correlation among  remaining variables are
#' within specified thresholds.
#' @param ldm LD correlation matrix
#' @param lcr lower correlation threshold
#' @param ucr upper correlation threshold
#' @param ... additional arguments (absorb and ignore)
#' @return index of variables to be retained.
ncv <- function(ldm, lcr=0, ucr=1, ...)
{
    if(lcr > 0 && ucr < 1)
    {
        r <- abs(ldm)
        r[upper.tri(r, TRUE)] <- lcr + (ucr - lcr) / 2
        b <- which(lcr <= r & r <= ucr, arr.ind=TRUE)
        apply(r, 2, function(.) all(lcr <= . & . <= ucr))
    }
    else
        rep(TRUE, ncol(ldm))
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
#' @param kgp.itr maximum iterations allowed (def=20)
#' @param drop return a vector if only one SNP is required (def=TRUE)
#' @param quiet report re-tries? (def=NO)
#' @param ... additional arguments
#'
#' Each time `kgp()` fails  to find a suitable set of  variants, it tries again
#' for at most `kgp.itr` times.
#' @export
kgp <- function(N, L=1,
                kgp.con=1, kgp.psd=NULL, kgp.ucr=NULL, kgp.maf=NULL,
                kgp.itr=20, drop=TRUE, quiet=TRUE, ...)
{
    psd <- kgp.psd %||% sqrt(.Machine$double.eps)
    ucr <- kgp.ucr %||% 0.99
    MAF <- kgp.maf %||% 0.05
    N17 <- nrow(C17)
    M17 <- ncol(C17)
    itr <- kgp.itr %||% 20

    ## L variants on demand, reserve 2 * L
    P <- min(L * 4, M17)
    while(itr)
    {
        i <- sample.int(N17, N, N > N17)          # N
        if(kgp.con)                               # P
            j <- seq(sample(M17 - P, 1) + 1, l=P) #
        else                                      #
            j <- sort(sample(M17, P))             #
        gmx <- C17[i, j, drop=FALSE]              # N x P
        gmx <- matrix(as.integer(gmx), N, P)      # int1 to int4
        gmx <- gmx[, maf(gmx) >= MAF, drop=FALSE] # min MAF
        gmx <- gmx[, asd(gmx) > 0, drop=FALSE]    # SD > 0
        if(NCOL(gmx) < L)
        {
            P <- min(P + L - NCOL(gmx), M17)
            if(!quiet)
                cat("P = ", P, "\n", sep="")
            next
        }
        
        ## when N > L, prevent colinearity
        if(N > L && psd > 0 && ucr < 1)
        {
            ldm <- cor(gmx)                      # P x P

            ## drop  SNP  to  enforce  non-linearity;  if there  are  less  than  L
            ## remaining, try again with a bigger reserve.
            gmx <- gmx[, ncv(ldm, ucr=ucr, ...), drop=FALSE]
            if(NCOL(gmx) < L)
            {
                P <- min(P + L - NCOL(gmx), M17)
                if(!quiet)
                    cat("P = ", P, "\n", sep="")
                next
            }
        }
            
        ## select L variants now
        j <- seq(sample(ncol(gmx) - L, 1) + 1, l=L)
        gmx <- gmx[, j, drop=drop]
            
        ## when N > L, prevent non-PD
        if(N > L && psd > 0 && ucr < 1)
        {
            ldm <- if(L > 1) cor(gmx) else 1
            ## if min(eigenvalue) < (threshold) * max(eigenvalue), try again
            egv <- eigen(ldm, TRUE, TRUE)$values
            if (egv[L] < psd * egv[1L])
            {
                if(!quiet)
                    cat("Non-PSD!\n")
                itr <- itr - 1
                next
            }
        }
        break
    }
    if(itr == 0)
    {
        msg <- gettextf("kgp() failed after %d tries with positive definite
            threshold kgp.psd=%.1e, relaxe kgp.psd to a positive value closer
            to 0 if N >= L, and a tiny negative if N < L.", kgp.itr, kgp.psd)
        stop(msg)
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
#' @param ... additional arguments (absorb and ignore)
#' @return N*L genotype matrix generated from two trail binomial
#' @export
bng <- function(N, L=1, MAF=NULL, drop=TRUE, std=FALSE, ...)
{
    MAF <- MAF %||% 0.05
    gmx <- stats::rbinom(N * L, 2L, MAF)
    if(std)
        gmx <- gmx / sqrt(2 * MAF * (1 - MAF))
    if(L > 1 || !drop)
        dim(gmx) <- c(N, L)
    gmx
}

#' Covernt real values to dosage genotype
#'
#' @param x matrix of read values, one variant per column;
#' @param m minor allele frequencies
#' @param ... additional arguments (absorb and ignore)
#' @return allele dosage in {0, 1, 2} with frequency {m^2, 2*m*(1-m), (1-m)^2}
#' @export
as.genotype <- function(x, m=NULL, ...)
{
    if (is.null(m))
    {
        m <- F17[sample(col(C17) - ncol(x), 1) + seq(ncol(x))]
    }
    for(j in seq(ncol(x)))
    {
        m2 <- m[j] * m[j]           # HWE freq of 0 allele
        mn <- 2 * m[j] * (1 - m[j]) # HWE freq of 1 allele
        qt <- stats::quantile(x[, j], c(m2, m2 + mn))
        x[, j] <- 2 - (x[, j] > qt[1]) - (x[, j] > qt[2])
    }
    x
}
