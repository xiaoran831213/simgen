
#' Multivariate negative square root of a PD matrix
#'
#' Decorrelate assocation test statistics between multiple phenotypes and g-variants.
#' 
#' @param C M x M matrix of correlation among M variants adjusted for covariants.
#' @param D Q x Q matrix of correlation among Q multivariate.
mns <- function(C, D=NULL, eps=NULL, ...)
{
    eps <- eps %||%  sqrt(.Machine$double.eps)
    D <- D %||% 1

    C <- eigen(C, TRUE)
    D <- eigen(D, TRUE)
    
    d <- kronecker(D$values, C$values)
    u <- kronecker(D$vectors, C$vectors)
    dim(d) <- NULL
    
    ## positive eigen values
    . <- d > max(d) * eps
    if(!all(.))
    {
        d <- d[  .]
        u <- u[, .]
    }
    L <- length(d)              # effective number of eigen
    
    ## square root
    d <- sqrt(1/d)
    H <- u %*% (d * t(u))       # U diag(d) U'
    H <- 0.5 * (H + t(H))
    
    list(H=H, L=L)
}

#' Multivariate Decorrelation by Orthogonal Transformation
#'
#' Decorrelate assocation test statistics between multiple phenotypes and g-variants.
#' 
#' @param Z Q x M matrix of z-scores between Q traits and M variants.
#' @param C M x M matrix of correlation among M variants adjusted for covariants.
#' @param D Q x Q matrix of correlation among Q traits.
mdt <- function(Z, C, D=NULL, ...)
{
    ret <- mns(C, D, ...)
    H <- ret$H
    L <- ret$L
    dim(Z) <- NULL
    X <- H %*% Z

    P <- 1 - pchisq(sum(X^2), L)
    c(list(X=X, P=P), ret)
}
