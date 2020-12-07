tmp1 <- function()
{
    sgm <- matrix(0, 7, 7)
    sgm[lower.tri(sgm)] <- rbeta(21, .2, 1)
    sgm <- sgm + t(sgm)
    diag(sgm) <- 1
    sgm <- as.matrix(nearPD(sgm, corr=TRUE)$mat)
    dat <- MASS::mvrnorm(100, rep(0, 7), sgm)

    g <- dat[, 1:5]
    h <- dat[, 6:7]

    cg1 <- cov2cor(solve(solve(cor(dat))[1:5, 1:5]))
    cg2 <- cor(lm(g ~ h)$resid)

    print(all.equal(cg1, cg2))
}


tmp2 <- function(N=100)
{
    ## high LD genotype, 5 SNPs
    LD <- matrix(0, 5, 5)
    LD[lower.tri(LD)] <- rbeta(10, 4, 1)
    LD <- LD + t(LD)
    diag(LD) <- 1
    LD <- as.matrix(nearPD(LD, corr=TRUE)$mat)
    g0 <- MASS::mvrnorm(N, rep(0, 5), LD)

    ## outcome, 2 dimensional correlated
    CY <- matrix(c(1, .5, .5, 1), 2, 2)
    y0 <- MASS::mvrnorm(N, rep(0, 2), CY)

    ## both SNPs and outcomes associate with covariates
    x1 <- y0 %*% matrix(rnorm(2 * 3), 2, 3)
    x2 <- g0 %*% matrix(rnorm(5 * 3), 5, 3)
    x0 <- x1 + x2 + matrix(rnorm(N * 3), N, 3) * 1

    ## residuals
    ry <- lm(y0 ~ x0)$resid
    rg <- lm(g0 ~ x0)$resid

    ## collider effect
    print(gwa_lm(ry ~ {rg})$pvl)
    print(gwa_lm(ry ~ {g0})$pvl)

    ## original effect
    print(gwa_lm(y0 ~ {g0})$pvl)
}
