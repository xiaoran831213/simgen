mtq <- function(Z, C, D=NULL, tol.egv=NULL, ...)
{
    tol.egv <- tol.egv %||%  sqrt(.Machine$double.eps)
    D <- D %||% 1
    dim(Z) <- NULL

    C <- eigen(C, TRUE, TRUE)$values
    D <- eigen(D, TRUE, TRUE)$values
    
    d <- kronecker(D, C)
    dim(d) <- NULL

    ## . <- d > max(d) * tol.egv
    ## if(!all(.))
    ##     d <- d[  .]
    L <- length(d)              # effective number of eigen

    ## P <- imhof(sum(Z^2), d, delta=rep(0, L))$Qq
    P <- davies(sum(Z^2), lambda = d)$Qq
    list(d=d, L=L, P=P)
}

pow <- function(rpt)
{
    rpt <- subset(rpt, se=-itr)
    grp <- subset(rpt, se=-c(pvl, egv, mcr))
    rpt <- by(rpt, grp, function(g)
    {
        cfg <- subset(g, se=-c(pvl, egv))[1, ]
        pow <- with(g, mean(pvl <= 0.05))
        egv <- with(g, mean(egv))
        mcr <- with(g, mean(mcr))
        cbind(cfg, pow=pow, egv=egv, rep=nrow(g))
    })
    rpt <- do.call(rbind, rpt)
    rpt
}
