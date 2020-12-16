## nsp <- function(C, D=NULL, tol.egv=NULL, ...)
## {
##     tol.egv <- tol.egv %||%  sqrt(.Machine$double.eps)
##     D <- D %||% 1

##     C <- eigen(C, TRUE)
##     D <- eigen(D, TRUE)

##     d1 <- D$values
##     d2 <- C$values
##     i1 <- d1 > max(d1) * tol.egv
##     i2 <- d2 > max(d2) * tol.egv
##     d1 <- d1[i1]
##     d2 <- d2[i2]
##     u1 <- D$vectors[, i1]
##     u2 <- C$vectors[, i2]
##     d <- kronecker(d1, d2)
##     u <- kronecker(u1, u2)
    
##     ## d <- kronecker(D$values, C$values)
##     ## u <- kronecker(D$vectors, C$vectors)
    
##     ## positive eigen values
##     ## . <- d > max(d) * tol.egv
##     ## if(!all(.))
##     ## {
##     ##     d <- d[  .]
##     ##     u <- u[, .]
##     ## }
##     L <- length(d)              # effective number of eigen
    
##     ## square root
##     dim(d) <- NULL
##     d <- sqrt(1/d)
##     H <- u %*% (d * t(u))       # U diag(d) U'
##     H <- 0.5 * (H + t(H))
    
##     list(H=H, L=L)
## }

tmp1 <- function()
{
    times = 5e3
    N <- 5e2
    M <- 3
    L <- 3

    pvl <- list()
    dfs <- list()
    for(i in seq(times))
    {
        print(i)
        ## G <- mvn(N, 0, sim_cor(L, 7, 1, tol.psd=1e-9))
        ## G <- mvn(N, 0, sim_cor(L, 1, 1, tol.psd=1e-9))
        Y <- std(matrix(rnorm(N * 2), N, 2))

        G <- get_gmx(N, L, psd=1e-8, ucr=.99, maf=.05)
        G <- fac(G)
        ## G <- matrix(rnorm(N * L), N, L)
        ## LHS <- cbind(Y^2, Y[, 1] * Y[, 2])
        LHS <- Y
        ## LHS <- ORQ(Y[, 1] * Y[, 2])
        ## LHS <- BCX(Y^2)
        mdl <- LHS ~ {G}
        gwa <- gwa_lm(mdl)
        C <- gwa_pld(mdl)
        D <- cor(gwa$rsp)
        
        . <- nsp(C, D, tol.egv=1e-8)
        H <- .$H
        z <- H %*% c(gwa$zsc)

        pvl[[i]] <- pchisq(sum(z^2), .$L, lower.tail = FALSE)
        dfs[[i]] <- .$L
    }
    pvl <- unlist(pvl)
    mean(pvl < 0.05)
}
