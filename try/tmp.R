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

tmp2 <- function(a=0, b=1, e=1, N=1e2, times=5e2)
{
    pvl <- replicate(times,
    {
        x <- rnorm(N)
        u <- rnorm(N)

        y1 <- rnorm(N, a * x + b * x * u) + e * rnorm(N)
        ## y1 <- rnorm(N, a * x, exp(b * x)) + e * rnorm(N)
        
        m1 <- lm(y1 ~ x)
        p1 <- summary(m1)$coef[2, 4]

        y2 <- resid(m1)^2
        m2 <- lm(y2 ~ x)
        p2 <- summary(m2)$coef[2, 4]

        c(p1=p1, p2=p2)
    })

    rowMeans(pvl < 0.05)
}


tmp3 <- function(a=0, b=1, e=1, times=5e2, N=1e2)
{
    ## w <- matrix(rnorm(N * N), N, N)
    ## w <- t(svd(w)$u)
    pvl <- replicate(times,
    {
        x <- rnorm(N)
        u <- rnorm(N)

        y1 <- rnorm(N, a * x + b * x * u) + e * rnorm(N)
        y1 <- bin(y1)

        m1 <- lm(y1 ~ x)
        h1 <- predict(m1, type='res')
        p1 <- summary(m1)$coef[2, 4]

        y2 <- log(resid(m1)^2)
        m2 <- lm(y2 ~ x)
        p2 <- summary(m2)$coef[2, 4]

        h2 <- h1 * (1 - h1)
        m3 <- lm(y2 ~ h2)
        p3 <- summary(m3)$coef[2, 4]
        c(p1=p1, p2=p2, p3=p3)
    })

    rowMeans(pvl < 0.05)
}

tmp4 <- function(a=0, b=1, e=1, times=5e2, N=1e2)
{
    pvl <- replicate(times,
    {
        x <- rnorm(N)
        u <- rnorm(N)

        y1 <- rnorm(N, a * x, exp(b * x)) + e * rnorm(N)
        
        m1 <- lm(y1 ~ x)
        p1 <- summary(m1)$coef[2, 4]
        
        y2 <- resid(m1)^2
        m2 <- lm(y2 ~ x)
        p2 <- summary(m2)$coef[2, 4]

        m3 <- try(glm(y2 ~ x, Gamma("log")))
        if(!inherits(m3, "try-error"))
        {
            p3 <- summary(m3)$coef[2, 4]
        }
        else
        {
            p3 <- NA
        }

        c(p1=p1, p2=p2, p3=p3)
    })

    rowMeans(pvl < 0.05, na.rm=TRUE)
}
