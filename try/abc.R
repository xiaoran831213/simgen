
abc <- function()
{
    ## Can we adjust for linear effect of x on y?

    ## Case 1: effect of x on y is linear
    set.seed(12345)

    x <- rep(seq(lv), each=1e3) # x as real number
    f <- as.factor(x)           # x as categories
    n <- length(x)
    Sd <- Mu <- rep(1, n <- length(x))

    ## for different variances between levels of x (H_A):
    ## for(i in 1:lv) Sd[which(x == i)] <- runif(1, 1, 2)
    y1 <- x + rnorm(n, 0, Sd) + rnorm(n, 0, 1)
    yc <- unlist(lapply(split(y1, f), scale))
    sapply(split(yc, f), mean) # group mean is 0
    sapply(split(yc, x), sd)   # group sd is 1
    y2 <- yc^2

    anova(lm(y2 ~ factor(x) + y1), lm(y2 ~ y1))[2, 6]  # factors capture y2, LR
    summary(lm(y2 ~ factor(x) + y1))$coef              # factors capture y2, Wald

    anova(lm(y2 ~ x + I(x^2) + y1), lm(y2 ~ y1))[2, 6] # real x and x^2 capture y2, LR
    summary(lm(y2 ~ x + I(x^2) + y1))$coef             # read x and x^2 capture y2, Wald
    
    anova(lm(y2 ~ x + y1), lm(y2 ~ y1))[2, 6] # real x do not capture y2, LR test
    summary(lm(y2 ~ x + y1))$coefficients     # real x do not capture y2, Wald test

    ## Case 2: x acts on y as a factor
    for(i in 1:lv)
        Mu[which(x == i)] <- runif(1, 0.5, 1)
    ## for different variances between levels of x (H_A):
    ## for(i in 1:lv) Sd[which(x == i)] <- runif(1, 1, 2)
    y2 <- ( y <- Mu + rnorm(n, sd=Sd) + rnorm(n) )^2
    for(i in 1:lv) cat(mean(y[which(x==i)]), var(y[which(x==i)]), "\n")
    ff <- lm(y2 ~ factor(x) + y)
    fr <- lm(y2 ~ y)
    anova(fr, ff)$"Pr(>F)"[2]
    summary(lm(y2 ~ x + y))$coefficients[2,4]

    ## Case 3: z is confounder, treating x as factor works
    library(gtools)
    N=1e5; s1=50; s2=10
    latent1=rnorm(N,sd=s1)
    z = latent1 + rnorm(N,sd=s2)
    x2 = latent1 + rnorm(N,sd=s2)
    y = z + rnorm(N,sd=s2)
    summary(lm(y ~ x2))$coefficients[2,4]
    summary(lm(y ~ x2 + z))$coefficients[2,4]
    fx <- quantcut(x2, q=7) # make a factor variable
    ff <- lm(y ~ fx + z)
    fr <- lm(y ~ z)
    anova(fr, ff)$"Pr(>F)"[2]
}

cs1 <- function()
{
    ## Can we adjust for linear effect of x on y?
    ## Case 1: effect of x on y is linear
    set.seed(12345)

    x1 <- rep(seq(lv), each=1e3) # x as real number
    f1 <- as.factor(x)           # x as categories
    x2 <- x1^2
    n <- length(x)
    Sd <- Mu <- rep(1, n <- length(x))

    ## for different variances between levels of x (H_A):
    ## for(i in 1:lv) Sd[which(x == i)] <- runif(1, 1, 2)
    y1 <- x1 + rnorm(n, 0, Sd) + rnorm(n, 0, 1) # x is not a vQTL on y1

    yc <- unlist(lapply(split(y1, x1), scale)) # center y by x-level
    sapply(split(yc, x1), mean)                # group mean is 0
    sapply(split(yc, x1), sd)                  # group sd is 1

    y2 <- yc^2

    anova(lm(y2 ~ factor(x) + y1), lm(y2 ~ y1))[2, 6]
    summary(lm(y2 ~ factor(x) + y1))$coef            

    anova(lm(y2 ~ x + I(x^2) + y1), lm(y2 ~ y1))[2, 6]
    summary(lm(y2 ~ x + I(x^2) + y1))$coef            
    
    anova(lm(y2 ~ x + y1), lm(y2 ~ y1))[2, 6]
    summary(lm(y2 ~ x + y1))$coefficients
}

