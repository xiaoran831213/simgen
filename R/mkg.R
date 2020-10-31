## make dosage genotype

#' Discretizing dosage by Hardy-Weinberg Equilibrium
#'
#' @param x matrix of allele dosages in continuous scale
#' @param m vector of allele frequencies
#' @return allele dosage in  {0, 1, 2}  with frequency {m^2, 2*m*(1-m), (1-m)^2}
#' @noRd
hwe <- function(x, m=NULL, ...)
{
    if (is.null(m))
        m <- runif(ncol(x), 0.05, 0.45)
    for(j in seq(ncol(x)))
    {
        m2 <- m[j] * m[j]           # HWE freq of 0 allele
        mn <- 2 * m[j] * (1 - m[j]) # HWE freq of 1 allele
        qt <- quantile(x[, j], c(m2, m2 + mn))
        x[, j] <- 2 - (x[, j] > qt[1]) - (x[, j] > qt[2])
    }
    x
}

#' Get a correlation matrix
#' @noRd
cmx <- function(L, alpha=0.9, beta=alpha, ...)
{
    ## correlation
    R <- matrix(0, L, L)
    R[lower.tri(R)] <- sort(2 * rbeta(L * (L - 1) / 2, alpha, beta) - 1)
    R <- R + t(R)
    diag(R) <- 1
    R <- fpd(R) # force PD
    cov2cor(R)
}

#' Make the right side
#'
#' @param N sample size
#' @param L number of g-variants
#' @param M number of covariates
#' @param P number of g-variants that is effective
#' @param e size of white noise
#' @return
#' a list,
#' \itemize{
#' \item{G}{N * L matrix of genotype.}
#' \item{X}{N * M matrix of covariates.}
#' \item{C}{D * D matrix of correlation, where D = M + L}
#' \item{R}{L * L matrix of genotype LD.}
#' \item{m}{L-vector of allele frequencies (or 0.5 * mean_dosage).}
#' }
sim <- function(nrep, expr, N=1e3, L=3, M=2, P=1, nsd=1, ...)
{
    cfg <- get.arg()
    cfg$expr <- NULL
    cfg$nrep <- NULL
    if(missing(expr))
    {
        expr <- quote(
        {
            val <- summary(Y)
            key <- names(val)
            val <- round(as.vector(val), 3)
            .d(mtd="sm5", key=key, val=val)
        })
    }
    ## correlation of all variables
    C <- cmx(L + M, ...)
    ## print(C)

    res <- list()
    for(i in seq(nrep))
    {
        ## L variants and M covariates
        X <- mvn(N, 0, C)          
        ## print(cor(X))

        ## variants discretized
        if(L > 0)
            X[, 1:L] <- hwe(X[, 1:L], ...)
        ## print(cor(X))
        G <- X[, 0 + seq(1, l=L), drop=FALSE]
        X <- X[, L + seq(1, l=M), drop=FALSE]
        
        ## Y = noise + genetics + environment
        Y <- rnorm(N, 0, nsd) # noise
        if(L > 0 && P > 0)         # genetics
            Y <- Y + std(G[, sample(L, P), drop=FALSE] %*% rnorm(P))
        if(M > 0) # enviornment
            Y <- Y + std(X %*% rnorm(M))
        Y <- drop(std(Y))

        val <- eval(substitute(expr))
        res[[i]] <- .d(itr=i, val)
    }
    res <- cbind(cfg, do.call(rbind, res))
    
    ## return
    res
}
