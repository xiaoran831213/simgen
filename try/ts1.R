ts1 <- function()
{
    ##
    library(CompQuadForm)
    library(MASS)
    library(Matrix)
    library(globaltest)

    ## Hardy-Weinberg Equilibrium quantiles for discretizing MVN to 0,1,2
    HWE.quantiles <- function(x, f) {
        q0 <- quantile(x,f)[1]
        q1 <- quantile(x,f)[2]
        i0 <- which(x < q0)
        i1 <- which(x >= q0 & x < q1)
        i2 <- which(x >= q1)
        x[i0] <- 0; x[i1] <- 1; x[i2] <- 2
        x
    }
    ## generate population of SNPs
    Pop.Size <- 1e5
    L <- 100
    R <- matrix(1, L, L)
    R[lower.tri(R)] <- sort(2*rbeta(L*(L-1)/2, 0.9, 0.9) - 1)
    R <- (R * lower.tri(R)) + t(R * lower.tri(R))
    diag(R) <- 1
    R <- as.matrix(nearPD(R, corr=TRUE, posd.tol = 1e-12, maxit = 1000)$mat)
    u <- 20*(2*(rbeta(L, 0.25, 0.25) - 0.5))
    R <- cov2cor(R + u %*% t(u))
    mm <- seq(from = 0.03, to = 0.97, length.out = L)
    d <- mvrnorm(n=Pop.Size, mu=mm, Sigma=R)
    ## get genotype frequencies and discretize MVN to 0,1,2
    ff <- list(); for(i in 1:L) ff[[i]] <- c(mm[i]^2, 2*mm[i]*(1-mm[i]) + mm[i]^2)
    for(i in 1:L) { d[,i] <- HWE.quantiles(d[,i], ff[[i]]) }
    min.freq <- 0.03 # minimum MAF
    
    Sz <- 1000 # sample size
    effects.SD <- 0.2
    n.eff  <- round(L/4) # number of true effects
    oo <- 100 # number of simulations

#### Sample effects from MVN
    rho <- -1/(n.eff-1) # minimum equicorr is bounded by -1/(n.eff-1)
    S <- (1-rho) * diag(1,n.eff,n.eff) + rho * matrix(1,n.eff,n.eff)
    u <- 10*(2*(rbeta(n.eff, 0.5, 0.5) - 0.5))
    R <- cov2cor(S + u %*% t(u)) # S
    Effects <- effects.SD * mvrnorm(n=oo, mu=rep(0, n.eff), Sigma=R)
    ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@10@"]]));##:ess-bp-end:##
    p.mr <- p.gt <- p.tq <- rep(NA, oo)
    ii <- 1
    for(ii in 1:oo) {
        eff <- rep(0, L)
        ix <- sample(seq(1:L), n.eff)
        eff[ix] <- Effects[ii,]
        repeat { # sample from the population
            ib <- sample(seq(1 : dim(d)[1]), Sz, replace=FALSE)
            gg <- as.matrix(d[ib,])
            f <- rep(NA, L <- dim(gg)[2])
            for(i in 1:L) f[i] <- (1 + mean(gg[,i]-1)) / 2 # 1-MAF, == sum(gg[,i]) / (dim(gg)[1]*2)
            if(min(f) >= min.freq) {
                Cr <- cor(gg)
                LD <- Cr[lower.tri(Cr)]
                if(max(abs(LD)) < 0.999) break
            }
        }
        X <- gg
        proxies <- seq(1:L)[-ix]
        fix.ix <- sample(proxies, min(10, length(proxies)))
        X.fixed <- gg[, fix.ix]
        X.random <- gg[, -fix.ix]
        LL <- dim(X.random)[2]
        rg <- t.val <- rep(0.99, LL)
        Y <- (X %*% eff) + rnorm(Sz, sd=1) # Outcome
        rg.fixed <- lm(Y ~ X.fixed)
        Y.res <- rg.fixed$residuals
        ## P-value for statistic Q
        ## first, convert X to data frame, because 'gt' needs colnames
        Xd <- data.frame(X)
        f.fixed <- as.formula(paste("Y ~ ", paste(names(Xd[fix.ix]), collapse= "+")))
        f.random <- as.formula(paste("Y ~ ", paste(names(Xd[-fix.ix]), collapse= "+")))
        p.gt[ii]  <- p.value( gt(f.fixed, f.random, data=Xd) )
        for(j in 1:LL) { # Marginal regressions
            sumlm <- summary(lm(Y.res ~ X.random[, j]))
            rg[j] <- sumlm$coefficients[2,4]
            t.val[j] <- sumlm$coefficients[2,3]
        }
        x <- summary(lm(Y.res ~ X.random))$fstatistic
        p.mr[ii] <- pf(x[1], x[2], x[3], lower.tail=F)
        Cx <- cor(X)
        Cr <- ginv(ginv(Cx)[-fix.ix, -fix.ix]) # use the Schur complement to adjust correlations
        eigva <- eigen(Cr, only.values = TRUE)$values
        davies(sum(t.val^2), lambda = eigva)$Qq
        p.tq[ii] <- davies(sum(t.val^2), lambda = eigva)$Qq
        if(ii %% 10 == 0) cat(".", ii, ".", sep="");
        if(ii %% 100 == 0) { cat("\r") }
    }
    cat("\n")
    par(mfrow=c(2,2));
    plot(-log(p.gt), -log(p.tq), xlim=c(0,30), ylim=c(0,30))
    plot(-log(p.gt), -log(p.mr), xlim=c(0,30), ylim=c(0,30))
    plot(p.gt, p.tq)
    plot(p.gt, p.mr)
    alpha <- 0.05
    cat("Tq, Q, and Multiple regression\n", ecdf(p.tq)(alpha), ecdf(p.gt)(alpha), 
        ecdf(p.mr)(alpha),"\n")
}
