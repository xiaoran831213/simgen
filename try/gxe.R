## test GxE models

## vQTL induced by interaction
gx1 <- function(a=0, b=1, e=1, times=5e2, N=1e2)
{
    pvl <- replicate(times,
    {
        x1 <- c(kgp(N, 1))
        u1 <- rnorm(N)
        x2 <- cbind(std(x1), std(std(x1)^2))
        
        y1 <- std(rnorm(N, a * x1 + b * u1 * x1, e))
        m1 <- lm(y1 ~ x1)
        p1 <- summary(m1)$coef[2, 4]

        y2 <- resid(m1)^2
        nl <- lm(y2 ~ 1)          # null model
        ma <- lm(y2 ~ x1)         # linear
        pa <- anova(ma, nl)[2, 6] # 
        mb <- lm(y2 ~ x2)         # quadratic
        pb <- anova(mb, nl)[2, 6] #

        ## DGLM
        mg <- try(glm(y2 ~ x1, Gamma("log")))
        pg <- if(inherits(mg, 'try-error')) NA else summary(mg)$coef[2, 4]
        
        c(p1=p1, pa=pa, pb=pb, pg=pg)
    })

    rowMeans(pvl < 0.05, na.rm=TRUE)
}

## vQTL incuded by parameterized variance
gx2 <- function(a=0, b=1, e=1, times=5e2, N=1e2)
{
    pvl <- replicate(times,
    {
        x1 <- c(kgp(N, 1))
        x1 <- std(x1, TRUE, FALSE)
        x2 <- std(cbind(x1, x1^2))

        y1 <- std(rnorm(N, a * x1, exp(b * x1) + e))
        m1 <- lm(y1 ~ x1)
        p1 <- summary(m1)$coef[2, 4]

        y2 <- resid(m1)^2
        nl <- lm(y2 ~ 1)          # null model
        ma <- lm(y2 ~ x1)         # linear
        pa <- anova(ma, nl)[2, 6] # 
        mb <- lm(y2 ~ x2)         # quadratic
        pb <- anova(mb, nl)[2, 6] #

        ## Gamma regression
        mg <- try(glm(y2 ~ x1, Gamma("log")))
        pg <- if(inherits(mg, 'try-error')) NA else summary(mg)$coef[2, 4]
        ## dt <- data.frame(y1=y1, x1=x1)
        ## mg <- try(dglm(y1 ~ x1, ~ x1, data=dt), TRUE)
        ## pg <- if(inherits(mg, 'try-error')) NA else summary(mg$dispersion.fit)$coef[2, 4]

        c(p1=p1, pa=pa, pb=pb, pg=pg)
    })

    rowMeans(pvl < 0.05, na.rm=TRUE)
}

## surrogate traits for a binary trait
bn1 <- function(a=0, b=1, e=1, times=5e2, N=1e2)
{
    ## w <- matrix(rnorm(N * N), N, N)
    ## w <- t(svd(w)$u)
    pvl <- replicate(times,
    {
        x1 <- std(rnorm(N))
        u1 <- std(rnorm(N))
        x2 <- std(x1^2)

        y1 <- bin(rnorm(N, a * x1 + b * u1 * x1, e))
        m1 <- glm(y1 ~ x1, 'binomial')
        p1 <- summary(m1)$coef[2, 4]

        y2 <- resid(m1)^2
        nl <- lm(y2 ~ 1)          # null
        ma <- lm(y2 ~ x1)         # linear
        pa <- anova(ma, nl)[2, 6] #
        mb <- lm(y2 ~ x1 + x2)    # quadratic
        pb <- anova(mb, nl)[2, 6] #

        ## Gamma regression
        mg <- try(dglm(y1 ~ x1, ~ x1, binomial("logit")), TRUE)
        pg <- if(inherits(mg, 'try-error')) NA else summary(mg$dispersion.fit)$coef[2, 4]

        c(p1=p1, pa=pa, pb=pb, pg=pg)
    })

    rowMeans(pvl < 0.05, na.rm=TRUE)
}

tmp4 <- function(a=0, b=1, e=1, times=5e2, N=1e2)
{
    ##An example using dglm.Pvalues(.) with 200 simulated observations
    set.seed(123)
    require(dglm)
    n = 200
    a <- rbinom(n, 1, 0.5) ##Covariate for additive genetic effects
    sex <- rbinom(n, 1, 0.5) ##Covariate for a non-genetic effect
    add.effect = 0
    sex.effect = 1
    res.var = exp( a*add.effect ) ##Residual variance
    y <- 10 + a*add.effect + sex*sex.effect + rnorm(n ,0 , sqrt(res.var))
    ##The additive genetic effect must be given first in the following formula
    d.fit <- dglm( formula=y~a+ sex, dformula=~a)
    P.values <- dglm.Pvalues( d.fit )
    print( P.values )
}

dglm.Pvalues <- function(dglm.fit)
{
    P.disp = anova.dglm(dglm.fit)$Adj.P[2]
    P.mean = summary(dglm.fit)$coef[2,4]
    list(P.mean=P.mean, P.disp=P.disp)
}
