# 
library(CompQuadForm)
library(MASS)
library(Matrix)
options(width = 180)

## set.seed(12345)

Sz <- 5e2 # sample size
oo <- 1e3 # number of simulations
L <- 5 # number of observed predictors
LL <- L + 0*L # number of observed + unobserved predictors
n.traits <- 5
n.conf <- max(L, n.traits)
# Generate correlation matrices for predictors and outcome
R <- matrix(1, LL, LL); R[lower.tri(R)] <- sort(2*rbeta(LL*(LL-1)/2, 0.25, 0.25) - 1)
R <- (R * lower.tri(R)) + t(R * lower.tri(R)); diag(R) <- 1
u <- 2*(2*(rbeta(LL, 0.25, 0.25) - 0.5))
R <- cov2cor(R + u %*% t(u))
RR <- as.matrix(nearPD(R, corr=TRUE, maxit=1000, posd.tol = 1e-8)$mat)
R <- RR[1:L, 1:L]
#
P <- matrix(1, n.traits, n.traits)
P[lower.tri(P)] <- sort(2*rbeta(n.traits*(n.traits-1)/2, 1, 1) - 1)
P <- (P * lower.tri(P)) + t(P * lower.tri(P)); diag(P) <- 1
u <- 2*(2*(rbeta(n.traits, 1, 1) - 0.5))
P <- cov2cor(P + u %*% t(u))
P <- as.matrix(nearPD(P, corr=TRUE, posd.tol = 1e-3)$mat)
# Correlation for confounders
S <- diag(n.conf)

p.tq <- p.dot <- p.tq.res <- p.dot.res <- rep(NA, oo)
min.eigv  <- 1e-14
#set.seed(1234)
ii <- 1
for(ii in 1:oo)
{
    Conf <- 1.25 * mvrnorm(Sz, mu=rep(0, n.conf), Sigma=S)
    X <- mvrnorm(n=Sz, mu=rep(0, LL), Sigma=RR)
    X[,1:L] <- X[,1:L] + Conf[,1:L]/n.conf * 2
    X <- scale(X)
    Y <- mvrnorm(n=Sz, mu=rep(0,n.traits), Sigma=P)
    Y <- Y + Conf[,1:n.traits]/n.conf * 2
    if(ii==1 && dim(Y)[2] < 9)
        print(cor(Y))

    X <- scale(X)
    Y <- scale(Y)
    Conf <- scale(Conf)
    nY <- dim(Y)[2]
    nX <- dim(X)[2]

    Cr <- ginv(ginv( cor(cbind(X,Conf)) )[1:nX, 1:nX])
    Y.res <- std(lm(Y~Conf)$resid)
    tt <- c(gwa_lm(Y ~ {X} + Conf)$zsc)
    tt.res <- c(gwa_lm(Y.res ~ {X})$zsc)

    e1 <- eigen(cor(Y), symmetric = TRUE)
    e2c <- eigen(cov2cor(Cr), symmetric = TRUE)
    eigva <- kronecker(e1$values, e2c$values)
    eivec <- kronecker(e1$vectors, e2c$vectors)
    ix <- seq(1:length(eigva))
    for(j in 1:length(eigva))
    {
        if(eigva[j] < min.eigv)
        {
            ix <- ix[-j]
        }
    }
    pc <- eivec[,ix] %*% diag(sqrt(1/eigva[ix])) %*% t(eivec[,ix])
    (p.tq[ii] <- davies(sum(tt^2), lambda = eigva)$Qq)
    (p.dot[ii] <- 1 - pchisq(sum( (tt %*% pc)^2 ), df = length(ix)))
    e1 <- eigen(cor(Y.res), symmetric = TRUE)
    e2 <- eigen(Cr, symmetric = TRUE)
    eigva <- kronecker(e1$values, e2$values)
    eivec <- kronecker(e1$vectors, e2$vectors)
    ix <- which(eigva >= min.eigv)
    ## ix <- seq(1:length(eigva))
    ## for(j in 1:length(eigva))
    ## {
    ##     if(eigva[j] < min.eigv)
    ##     {
    ##         ix <- ix[-j]
    ##     }
    ## }
    pc <- eivec[,ix] %*% diag(sqrt(1/eigva[ix])) %*% t(eivec[,ix])
    (p.tq.res[ii] <- davies(sum(tt.res^2), lambda = eigva)$Qq)
    (p.dot.res[ii] <- 1 - pchisq(sum( (tt.res %*% pc)^2 ), df = length(ix)))
    if(!(ii %% 10)) cat(".", ii, ".", sep=""); if(!(ii %% 100)) cat("\r")
}
cat("\n")
a <- 0.05
cat(ecdf(p.tq)(a),  ecdf(p.tq.res)(a), ecdf(p.dot)(a), ecdf(p.dot.res)(a), "\n")
mean(diag(cor(Y, Conf))); mean(diag(cor(X, Conf)))

par(mfrow=c(1,2))
plot(p.tq, p.tq.res)
plot(p.dot, p.dot.res)
