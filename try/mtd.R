mtq <- function(Z, C, D=NULL, eps=NULL, ...)
{
    eps <- eps %||%  sqrt(.Machine$double.eps)
    D <- D %||% 1
    dim(Z) <- NULL

    C <- eigen(C, TRUE, TRUE)$values
    D <- eigen(D, TRUE, TRUE)$values
    
    d <- kronecker(D, C)
    dim(d) <- NULL

    . <- d > max(d) * eps
    if(!all(.))
        d <- d[  .]
    L <- length(d)              # effective number of eigen

    P <- davies(sum(Z^2), lambda = d)$Qq
    list(d=d, L=L, P=P)
}
