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

        G <- kgp(N, L, psd=1e-8, ucr=.99, maf=.05)
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

library(microbenchmark)
tmp2 <- function(N=5e2, M=3, a=0, b=1, e=2)
{
    rhm <- ~G(5) + X(5) | LD(G, a=6, b=1) + CX(X, a=1, b=1) + .5 @ GX(G+X, a=1, b=1)
    lhm <- Y(M) ~ a @ G + b @ X 
    flood(sim_dsg(rhm, N)$dat)
    flood(sim_rsp(lhm)$dat)
    Y <- Y + matrix(rnorm(N * M), N, M) * e # outcomes
    Y <- std(Y, TRUE, FALSE)
    X <- std(X, TRUE, FALSE)
    G <- std(G, TRUE, FALSE)
    Z <- get_rsd(Y, X, int=FALSE)
    H <- get_rsd(G, X, int=FALSE)
    m1 <- lm(Y ~ G + X)
    m2 <- lm(Z ~ H)

    s1 <- summary(m1)[[1]]$coef
    s2 <- summary(m2)[[1]]$coef
    list(s1, s2)
}
