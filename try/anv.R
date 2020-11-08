#' Check bugs in anova test
#'
#' @param N  # of samples
#' @param L  # of variants
#' @param P  # of casual variants
#' @param M  # of covariates
#' @param nsd SD of noise
main <- function(N=1e2, L=5, M=5, P=2, nsd=5)
{
    ## simulate dosage data
    R <- cmx(L + M)   # correlation among all variables
    X <- mvn(N, 0, R) # data matrix of variables
    A <- X[, -(1:L)]  # covariate -- fix effects
    B <- X[, +(1:L)]  # g-variate -- rnd effects
    B <- hwe(B)       # discretize by HWE

    ## effective variants
    Y <- rnorm(N, 0, nsd)             # noise
    Y <- Y + std(A %*% rnorm(M))      # covariate, effect
    p <- sample(L, P)                 # g-variate, casual ones
    Y <- Y + std(B[, p] %*% rnorm(P)) # g-variate, effect
    y <- drop(std(Y))

    ## models
    yA_ <- lm(Y ~ A)     # Y ~ A
    yAB <- lm(Y ~ A + B) # Y ~ A + B
    r <- resid(yA_)      # Y - A
    rA_ <- lm(r ~ A)     # Y - A ~ A
    rAB <- lm(r ~ A + B) # Y - A ~ A + B

    ## compare residuals
    print(paste("residual of Y ~ A     eq  (Y - A) ~ A    :", all.equal(resid(yA_), resid(rA_))))
    print(paste("residual of Y ~ A + B eq  (Y - A) ~ A + B:", all.equal(resid(yAB), resid(rAB))))
    print("")

    ## compare models
    rbind(.d(LRT="yAB vs yA_", anova(yAB, yA_)[2, c(3, 5, 6)]),
          .d(LRT="rAB vs rA_", anova(rAB, rA_)[2, c(3, 5, 6)]))
}
