#' Temporary Test
#'
#' @param N  # of samples
#' @param L  # of variants
#' @param P  # of casual variants
#' @param M  # of covariates
#' @param nsd SD of noise
ts3 <- function(N=1e2, L=5, M=5, P=2, nsd=5)
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
    p.tq <- p.dot <- p.tq.res <- p.dot.res <- rep(NA, oo)
    p.tq2 <- p.dot2 <- p.tq2.res <- p.dot2.res <- rep(NA, oo)
    ##
    ##set.seed(1234)
    ii <- 1
    for(ii in 1:oo)
    {
        Conf <- mvn(ssz, 0, S)                    # confounders
        X <- mvn(ssz, 0, RR)                      # SNPs
        X[, 1:L] <- X[, 1:L] + Conf[, 1:L]/n.conf # confounders on SNPs
        X[, 1:L] <- hwe(X[, 1:L])                 # discretize SNPs
        X <- std(X)                               #
        A <- X[, (L+1):LL]                        # unobserbed
        B <-  X[,1:L]                             # SNPs, observed
        X.prod <- model.matrix(~B:A-1)            # 2-way interaction
        np <- dim(X.prod)[2]

        ## traits = white noise + confounder effect
        Y <- mvn(ssz, 0, P) + Conf[,1:n.traits]/n.conf
        if(1)
        {
            for(i in 1:n.traits)
            {
                eff.prod <- 0.05 * runif(np, -1, 1)
                effects <- 0.05 * runif(LL, -1, 1) ## effects <- rep(0, LL)    
                Y[,i] <- Y[,i] +
                    0.5*c(X %*% effects) +
                    0.5*c(X.prod %*% eff.prod) +
                    0.5*sqrt(c(X^2 %*% abs(effects)))
                ##Y[,i] <- Y[,i] + 0.5*c(abs(X) %*% abs(effects))
                ##Y[,i] <- Y[,i] + 0.2*c(X %*% effects) + sqrt(c(X^2 %*% abs(effects)))
            }
        }
        if(ii==1) print(cor(Y))

        ## regressions
        X <- X[,1:L] ## drop "unobserved" part (X == B)
        Conf <- cbind(Conf, Y)
        Y <- cbind(.p(Y, 2, 1), .p(Y, 2, 2))
        cor(Y)
        X <- std(X)
        Y <- std(Y)
        ## Swap Y and X?
        ##Y.tmp  <- Y; X.tmp <- X; Y <- X.tmp; X <- Y.tmp
        ##
        nY <- dim(Y)[2]
        nX <- dim(X)[2]

        ## correlation among SNP X, but complement of confounders 
        Cr <- ginv(ginv( cor(cbind(X,Conf)) )[1:nX, 1:nX])
        tt.res  <- tt  <- t.check <- rep(0, nY*nX)
        tt2.res <- tt2 <- rep(0, n.traits*nX)
        v <- u <- 1
        for(i in 1:nY)
        {
            Y.res <- lm(Y[,i] ~ Conf)$residuals
            for(j in 1:nX)
            {
                tt[v] <- summary(lm(Y[,i] ~ X[, j] + Conf))$coefficients[2,3]
                ## NOTE: tt.res[v] --> sqrt(Cr[j,j])*tt[v]
                tt.res[v] <- summary(lm(Y.res ~ X[, j]))$coefficients[2,3]
                t.check[v] <- sqrt(Cr[j,j])*tt[v]
                if(i <= n.traits) {
                    tt2[u] <- tt[v]
                    tt2.res[u] <- tt.res[v]
                    u <- u+1
                }
                v <- v+1
            }
        }
        print(all.equal(matrix(tt,     nX, nY), gwa_exe(Y ~ {X} + Conf)))
        print(all.equal(matrix(tt.res, nX, nY), gwa_exe(Y - Conf ~ {X})))

        e1 <- eigen(cor(Y), symmetric = TRUE)
        e2c <- eigen(cov2cor(Cr), symmetric = TRUE)
        eigva <- kronecker(e1$values, e2c$values)
        eivec <- kronecker(e1$vectors, e2c$vectors)
        sD <- diag(sqrt(1/eigva))
        pc <- eivec %*% sD %*% t(eivec)
        p.tq[ii] <- davies(sum(tt^2), lambda = eigva)$Qq
        p.dot[ii] <- 1 - pchisq(sum( (tt %*% pc)^2 ), df = nY*nX)
        e2 <- eigen(Cr, symmetric = TRUE)
        eigva <- kronecker(e1$values, e2$values)
        eivec <- kronecker(e1$vectors, e2$vectors)
        sD <- diag(sqrt(1/eigva))
        pc <- eivec %*% sD %*% t(eivec)
        p.tq.res[ii] <- davies(sum(tt.res^2), lambda = eigva)$Qq
        p.dot.res[ii] <- 1 - pchisq(sum( (tt.res %*% pc)^2 ), df = nY*nX)
        ## analysis without Yi*Yj
        nY <- n.traits ## use only Y^2, drop Yi*Yj
        e1 <- eigen(cor(Y[,1:nY]), symmetric = TRUE)
        eigva <- kronecker(e1$values, e2c$values)
        eivec <- kronecker(e1$vectors, e2c$vectors)
        sD <- diag(sqrt(1/eigva))
        pc <- eivec %*% sD %*% t(eivec)
        p.tq2[ii] <- davies(sum(tt2^2), lambda = eigva)$Qq
        p.dot2[ii] <- 1 - pchisq(sum( (tt2 %*% pc)^2 ), df = nY*nX)
        eigva <- kronecker(e1$values, e2$values)
        eivec <- kronecker(e1$vectors, e2$vectors)
        sD <- diag(sqrt(1/eigva))
        pc <- eivec %*% sD %*% t(eivec)
        p.tq2.res[ii] <- davies(sum(tt2.res^2), lambda = eigva)$Qq
        p.dot2.res[ii] <- 1 - pchisq(sum( (tt2.res %*% pc)^2 ), df = nY*nX)
        if(!(ii %% 10)) cat(".", ii, ".", sep=""); if(!(ii %% 100)) cat("\r")
    }
    cat("\n")
    cat(ecdf(p.tq)(0.05),  ecdf(p.tq.res)(0.05), ecdf(p.dot)(0.05), ecdf(p.dot.res)(0.05), "\n")
    cat(ecdf(p.tq2)(0.05), ecdf(p.tq2.res)(0.05), ecdf(p.dot2)(0.05), ecdf(p.dot2.res)(0.05), "\n")

    (((ssz-1) / ssz) * tt.res) / t.check
    (((ssz-2) / ssz) * tt.res) / t.check
    (((ssz-3) / ssz) * tt.res) / t.check
    (((ssz-4) / ssz) * tt.res) / t.check
    1 / ( (((ssz-5) / ssz) * tt.res) / t.check )
    1 / ( (((ssz-6) / ssz) * tt.res) / t.check )


    ##x1 <- p.tq.res
    ##x2 <- p.tq
    ##x3 <- p.dot.res
    ##x4 <- p.dot

    ##par(mfrow=c(2,2))
    ##plot(x1 , p.tq.res)
    ##plot(x2 , p.tq)
    ##plot(x3 , p.dot.res)
    ##plot(x4 , p.dot)

    ##par(mfrow=c(2,2))
    ##plot(p.tq2.res , p.tq.res)
    ##plot(p.tq2 , p.tq)
    ##plot(p.dot2.res , p.dot.res)
    ##plot(p.dot2 , p.dot)
}
