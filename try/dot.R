dot <- function(Z, C, tol.cor=NULL, tol.egv=NULL, ...)
{
    if(is.null(tol.cor))
        tol.cor <- sqrt(.Machine$double.eps)
    if(is.null(tol.egv))
        tol.egv <- sqrt(.Machine$double.eps)

    ## trim collinear variants
    m <- dvt(C, tol.cor)
    M <- sum(m)                      # effective number of variants
    C <- C[m, m]
    Z <- Z[m]

    ## get orthogonal transformation
    d <- nsp(C, eps=tol.egv, ...)
    H <- d$H                         # orthogonal transformation
    L <- d$L                         # effective number of eigenvalues
    
    X <- H %*% Z                     # decorrelated statistics

    list(H=H, X=X, L=L, M=M)
}

dot_chisq <- function(Z, C, ...)
{
    ret <- dot(Z, C, ...)        # decorrelate
    L <- ret$L                   # effective number of eigenvalues

    Y <- sum(ret$X^2)            # sum of squares
    P <- 1 - pchisq(Y, L)        # a single p-value
    c(list(P=P, Y=Y), ret)
}
