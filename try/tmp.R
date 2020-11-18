#' Temporary Test
#'
#' @param N  # of samples
#' @param L  # of variants
#' @param P  # of casual variants
#' @param M  # of covariates
#' @param nsd SD of noise
main <- function(N=1e2, L=5, M=5, P=2, nsd=5)
{
    library(CompQuadForm)
    library(MASS)
    library(Matrix)

    ## Include.Cov=1 currently works with Tq, but not DOT!!!!
    Include.Cov <- 1

    ## set.seed(12345)
    HWE.quantiles <- function(x, f) {
        q0 <- quantile(x,f)[1]; q1 <- quantile(x,f)[2]
        i0 <- which(x < q0)
        i1 <- which(x >= q0 & x < q1)
        i2 <- which(x >= q1)
        x[i0] <- 0; x[i1] <- 1; x[i2] <- 2; x
    }

    Sz <- 1000 # sample size
    oo <- 100 # number of simulations
    Non.linear.effects  <- 1
    L <- 4 # number of observed predictors
    LL <- L + 25*L # number of observed + unobserved predictors
    n.traits <- 5
    n.conf <- max(L, n.traits)

    ## Generate correlation matrices for predictors and outcome
    R <- matrix(1, LL, LL); R[lower.tri(R)] <- sort(2*rbeta(LL*(LL-1)/2, 0.25, 0.25) - 1)
    R <- (R * lower.tri(R)) + t(R * lower.tri(R)); diag(R) <- 1
    u <- 2*(2*(rbeta(LL, 0.25, 0.25) - 0.5))
    R <- cov2cor(R + u %*% t(u))
    RR <- as.matrix(nearPD(R, corr=TRUE, maxit=1000, posd.tol = 1e-5)$mat)
    R <- RR[1:L, 1:L]
                                        #
    P <- matrix(1, n.traits, n.traits)
    P[lower.tri(P)] <- sort(2*rbeta(n.traits*(n.traits-1)/2, 1, 1) - 1)
    P <- (P * lower.tri(P)) + t(P * lower.tri(P)); diag(P) <- 1
    u <- 2*(2*(rbeta(n.traits, 1, 1) - 0.5))
    P <- cov2cor(P + u %*% t(u))
    P <- as.matrix(nearPD(P, corr=TRUE, posd.tol = 1e-5)$mat)
                                        # Correlation for confounders
    S <- matrix(1, n.conf, n.conf)
                                        #S[lower.tri(S)] <- 0
    S[lower.tri(S)] <- sort(2*rbeta(n.conf*(n.conf-1)/2, 1, 1) - 1)
    S <- (S * lower.tri(S)) + t(S * lower.tri(S)); diag(S) <- 1
    u <- 2*(2*(rbeta(n.conf, 1, 1) - 0.5))
    S <- cov2cor(S + u %*% t(u))
    S <- as.matrix(nearPD(S, corr=TRUE, posd.tol = 1e-3)$mat)
                                        #
    p.tq <- p.dot <- p.tq.res <- p.dot.res <- rep(NA, oo)
    p.tq.var <- p.dot.var <- p.tq.var.res <- p.dot.var.res <- rep(NA, oo)
    p.tq.cov <- p.dot.cov <- p.tq.cov.res <- p.dot.cov.res <- rep(NA, oo)
                                        #set.seed(1234)
    ii <- 1
    for(ii in 1:oo)
    {
        Conf <- mvrnorm(n=Sz, mu=rep(0, n.conf), Sigma=S) # confounders
        X <- mvrnorm(n=Sz, mu=rep(0, LL), Sigma=RR)
        X[,1:L] <- X[,1:L] + Conf[,1:L]/n.conf
        MAF <- runif(L, 0.05, 0.95)
        ff <- list(); for(i in 1:L) ff[[i]] <- c(MAF[i]^2, 2*MAF[i]*(1-MAF[i]) + MAF[i]^2)
        for(i in 1:L) { X[,i] <- HWE.quantiles(X[,i], ff[[i]]) }
        X <- scale(X)
        A <- X[, (L+1):LL]
        B <-  X[,1:L]
        X.prod <- A[, rep(seq(ncol(A)), each=ncol(B))] * B[, rep(seq(ncol(B)), ncol(A))] # kronecker on rows of A,B
        np <- dim(X.prod)[2]
        Y <- mvrnorm(n=Sz, mu=rep(0,n.traits), Sigma=P) + Conf[,1:n.traits]/n.conf
        if(0) { # <-- HA or H0?
            for(i in 1:n.traits) {
                eff.prod <- 0.1*runif(np, -1, 1)
                effects <- 0.1*runif(LL, -1, 1)
                if(Non.linear.effects) {
                    Y[,i] <- Y[,i] + 0.1*c(X %*% effects) + 0.25*c(X.prod %*% eff.prod) + 0.5*sqrt(c(X^2 %*% abs(effects)))
                } else {
                    Y[,i] <- Y[,i] + c(X %*% effects)
                }
            }
        }
        if(ii==1 && dim(Y)[2] < 9) print(cor(Y))
        X <- X[,1:L] # drop "unobserved" part
        X <- scale(X)
        Y <- scale(Y)
        nX <- dim(X)[2]
        Cr <- ginv(ginv( cov(cbind(X,Conf)) )[1:nX, 1:nX])
        Z <- NULL;

        for(i in 1 : dim(Y)[2]) {
            for(j in 1:nX) {
                Z <- cbind(Z, lm(Y[,i] ~ X[,j])$residuals)
            }
        }
        Z <- scale(Z)
        nZ <- dim(Z)[2]
        nY <- n.traits
        Z <- cbind(Z, Z^2)
        if(Include.Cov) {
            for(i in 1 : (n.traits-1)) {
                for(j in (i+1) : n.traits) {
                    nY <- nY + 1
                    for(k in 1:nX) {
                        Z <- cbind(Z, Z[,i+k-1]*Z[,j+k-1])
                    }
                }
            }
        }
        Z <- Z[, -c(1 : nZ)]
        Z <- scale(Z)
        ld <- kronecker(matrix(1, nrow = dim(Z)[2]/dim(Cr)[2], ncol = dim(Z)[2]/dim(Cr)[2]), Cr)
        Zcov <- cor(Z) * ld
        if(ii==1) {
            Zcov.all <- Zcov; ttm <- matrix(nrow=oo, ncol=nY*nX)
        } else {
            Zcov.all <- Zcov.all + Zcov
        }
        tt <- rep(0, nY*nX)
        v <- 1
        for(i in 1:nY) {
            for(j in 1:nX) {
                tt[v] <- summary(lm(Z[,v] ~ X[,j] + Conf))$coefficients[2,3]
                v <- v+1
            }
        }
        ttm[ii,] <- tt
        ee <- eigen(cov2cor(Zcov), symmetric = TRUE)
        eigva <- ee$values
        eivec <- ee$vectors
        (p.tq[ii] <- davies(sum(tt^2), lambda = eigva)$Qq)
        Le <- length(eigva)
        for(j in 2:length(eigva)) { if(eigva[j] < 5e-2) { Le <- j-1; break }}
        eigva <- eigva[1 : Le]
        eivec <- eivec[, 1 : Le]
        sD <- diag(sqrt(1/eigva))
        pc <- eivec %*% sD %*% t(eivec)
        (p.dot[ii] <- 1 - pchisq(sum( (tt %*% pc)^2 ), df = Le))
        if(!(ii %% 10)) cat(".", ii, ".", sep=""); if(!(ii %% 100)) cat("\r")
    }
    cat("\n")
    Zcor <- cov2cor(Zcov.all/oo)
    a <- 0.05
    cat(ecdf(p.tq)(a),  ecdf(p.dot)(a), "\n")

    plot(p.tq, p.dot)

    i <- min(12, dim(ttm)[2])
    cor(ttm)[1:i,1:i]; Zcor[1:i,1:i]
}
