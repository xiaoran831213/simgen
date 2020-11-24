mtq <- function(Z, C, D=NULL, tol.egv=NULL, ...)
{
    tol.egv <- tol.egv %||%  sqrt(.Machine$double.eps)
    D <- D %||% 1
    dim(Z) <- NULL

    C <- eigen(C, TRUE, TRUE)$values
    D <- eigen(D, TRUE, TRUE)$values
    
    d <- kronecker(D, C)
    dim(d) <- NULL

    . <- d > max(d) * tol.egv
    if(!all(.))
        d <- d[  .]
    L <- length(d)              # effective number of eigen

    Q <- imhof(sum(Z^2), d, delta=rep(0, L))$Qq
    ## P <- davies(sum(Z^2), lambda = d)$Qq
    list(d=d, L=L, P=Q)
}
