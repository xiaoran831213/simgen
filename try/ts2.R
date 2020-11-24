#' Temporary Test
#'
#' @param N  # of samples
#' @param L  # of variants
#' @param P  # of casual variants
#' @param M  # of covariates
#' @param nsd SD of noise
ts2 <- function(N=1e2, L=5, M=5, P=2, nsd=5)
{
    ## v2: take Y^2, Yi*Yj, rather than X^2, Xi*Xj
    library(CompQuadForm)
    library(MASS)
    library(Matrix)
    ## library(globaltest)

    set.seed(1234)

    ssz <- 1000 ## sample size
    oo <- 100 ## number of simulations
    L <- 5 ## number of observed predictors
    LL <- L + 20*L ## number of observed + unobserved predictors
    n.traits <- 4
    n.conf <- max(L, n.traits)
    ## Generate correlation matrices for predictors and outcome
    RR <- cmx(LL, rho=log(LL)) # all predictors
    R <- RR[1:L, 1:L]          # observed SNPs
    P <- cmx(n.traits)         # traits
    S <- cmx(n.conf)           # confounders
    ##
    p.tq <- p.dt <- p.tq.res <- p.dt.res <- rep(NA, oo)
    p.tq2 <- p.dt2 <- p.tq2.res <- p.dt2.res <- rep(NA, oo)
    ##
    ##set.seed(1234)
    ii <- 1
    for(ii in 1:oo)
    {
        Conf <- mvn(ssz, 0, S)                      # confounders
        X <- mvn(ssz, 0, RR)                        # SNPs
        X[, 1:L] <- X[, 1:L] + Conf[, 1:L] / n.conf # confounders on SNPs
        X[, 1:L] <- hwe(X[, 1:L])                   # discretize SNPs
        X <- std(X)                                 #
        A <- X[,(L+1):LL]                           # unobserbed
        B <-  X[,1:L]                               # SNPs, observed
        X.prod <- frm_mtx(~B:A)
        np <- dim(X.prod)[2]

        ## traits = white noise + confounder effect
        Y <- mvn(ssz, 0, P) + Conf[,1:n.traits]/n.conf
        if(0)
        {
            for(i in 1:n.traits)
            {
                eff.prod <- 0.05 * runif(np, -1, 1)
                effects <- 0.05 * runif(LL, -1, 1) ## effects <- rep(0, LL)    
                Y[,i] <- Y[,i] +
                    0.5*c(X %*% effects) +
                    0.5*c(X.prod %*% eff.prod) +
                    0.5*sqrt(c(X^2 %*% abs(effects)))
            }
        }
        if(ii==1) print(cor(Y))

        ## regressions
        X <- X[,1:L] ## drop "unobserved" part (X == B)
        Conf <- cbind(Conf, Y)
        Y <- cbind(.p(Y, 2, 1), .p(Y, 2, 2))
        X <- std(X)
        Y <- std(Y)

        ## cor among Y and X, the latter is controlled of Conf
        Cr <- gwa_cst(Y ~ Conf + {X})$gmx
        
        ## GWAS between nY and nX
        tt1 <- gwa_lm(Y ~ Conf + {X})  # original t tests
        rt1 <- gwa_lm(Y - Conf ~ {X})  # residual t tests
        tt2 <- tt1[, 1:n.traits]       # only Y_i^2, no Y_i Y_j
        rt2 <- rt1[, 1:n.traits]       # only Y_i^2, no Y_i Y_j
        tc <- sqrt(diag(Cr)) * tt1

        dim(tt1) <- NULL
        dim(rt1) <- NULL
        dim(tt2) <- NULL
        dim(rt2) <- NULL

        ## decorrelate test statistics
        p.tq[ii] <- mtq(tt1, cov2cor(Cr), cor(Y))$P
        p.dt[ii] <- mdt(tt1, cov2cor(Cr), cov(Y))$P
        
        p.tq.res[ii] <- mtq(rt1, Cr, cor(Y))$P
        p.dt.res[ii] <- mdt(rt1, Cr, cor(Y))$P

        ## analysis without Yi*Yj
        p.tq2[ii] <- mtq(tt2, cov2cor(Cr), cor(Y[, 1:n.traits]))$P
        p.dt2[ii] <- mdt(tt2, cov2cor(Cr), cor(Y[, 1:n.traits]))$P

        p.tq2.res[ii] <- mtq(rt2, Cr, cor(Y[, 1:n.traits]))$P
        p.dt2.res[ii] <- mdt(rt2, Cr, cor(Y[, 1:n.traits]))$P

        if(!(ii %% 10)) cat(".", ii, ".", sep=""); if(!(ii %% 100)) cat("\r")
    }
    cat("\n")
    cat(ecdf(p.tq)(0.05),  ecdf(p.tq.res)(0.05), ecdf(p.dt)(0.05), ecdf(p.dt.res)(0.05), "\n")
    cat(ecdf(p.tq2)(0.05), ecdf(p.tq2.res)(0.05), ecdf(p.dt2)(0.05), ecdf(p.dt2.res)(0.05), "\n")

    (((ssz-1) / ssz) * rt1) / tc
    (((ssz-2) / ssz) * rt1) / tc
    (((ssz-3) / ssz) * rt1) / tc
    (((ssz-4) / ssz) * rt1) / tc
    1 / ( (((ssz-5) / ssz) * rt1) / tc )
    1 / ( (((ssz-6) / ssz) * rt1) / tc )


    ##x1 <- p.tq.res
    ##x2 <- p.tq
    ##x3 <- p.dt.res
    ##x4 <- p.dt

    ##par(mfrow=c(2,2))
    ##plot(x1 , p.tq.res)
    ##plot(x2 , p.tq)
    ##plot(x3 , p.dt.res)
    ##plot(x4 , p.dt)

    ##par(mfrow=c(2,2))
    ##plot(p.tq2.res , p.tq.res)
    ##plot(p.tq2 , p.tq)
    ##plot(p.dt2.res , p.dt.res)
    ##plot(p.dt2 , p.dt)
}
