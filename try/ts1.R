#' Temporary Test
#'
#' @param N  # of samples
#' @param M  # of traits
#' @param nsd SD of noise
tst1 <- function(N=2e2, M=2, a=0, b=1, e=1, ydt=0, evt=1e-4, seed=NULL, times=1e3, ...)
{
    arg <- get.arg(skp=c("seed", "out", "times"))
    psd <- arg$psd %||% evt
    set.seed(seed)

    rpt <- list()
    rhs <- ~G(5) + X(5) | HL(G, a=1, b=1) + LC(X, a=1, b=1)
    for(ii in seq(times))
    {
        flood(sim_dsg(rhs, N, psd=psd)$dat)
        flood(sim_rsp(Y(M) ~ a @ G + b @ X)$dat)
        Y <- Y + matrix(rnorm(N * M, 0, e), N, M)
        ## residual of Y ~ X
        Z <- as.matrix(resid(lm(Y ~ X)))
        ## residual of G ~ X
        J <- as.matrix(resid(lm(G ~ X)))

        cov_J_ <- scp(cor(cbind(G, X)), "X")
        ## stopifnot(all.equal(cov_J_, cov(J))) # not true
        cor_J_ <- cov2cor(cov_J_)
        ## stopifnot(all.equal(cov_J_, cor(J))) # true
        
        ## de-correlate
        if(ydt)
        {
            Z <- pcs(Z)
            Y <- pcs(Y)
        }

        ## GWAS between nY and nX
        ## t0 <- gwa_lm(Y ~     {G}) # no covariate
        t0 <- gwa_lm(Y ~ X + {G}) # original t tests
        t1 <- gwa_lm(Z ~     {G}) # residual t tests
        t2 <- gwa_lm(Z ~     {J}) # double residual

        ## test statistics
        r <- list()
        r[['D0a']] <- mdt(t0$zsc, cor(J), cor(Y), tol.egv=evt) #
        r[['D1a']] <- mdt(t1$zsc, cov_J_, cor(Z), tol.egv=evt) # controlled
        r[['D2a']] <- mdt(t2$zsc, cor(J), cor(Z), tol.egv=evt) # controlled
        r[['D2b']] <- mdt(t2$zsc, cor_J_, cor(Z), tol.egv=evt) # controlled

        r[['T0a']] <- mdt(t0$zsc, cor(J), cor(Y), tol.egv=evt) # controlled
        r[['T1a']] <- mdt(t1$zsc, cov_J_, cor(Z), tol.egv=evt) # controlled
        r[['T2a']] <- mdt(t2$zsc, cor(J), cor(Z), tol.egv=evt) # controlled
        r[['T2b']] <- mdt(t2$zsc, cor_J_, cor(Z), tol.egv=evt) # controlled

        p <- sapply(r, `[[`, 'P')
        l <- sapply(r, `[[`, 'L')

        rpt[[ii]] <- .d(itr=ii, mtd=names(r), pvl=p, egv=l)
        cgy[[ii]] <- mean(abs(cor(G, Y)))
        cgx[[ii]] <- mean(abs(cor(G, X)))
        cxy[[ii]] <- mean(abs(cor(X, Y)))
        if(!(ii %% 10))
        {
            cat(".", ii, ".", sep=""); if(!(ii %% 100)) cat("\r")
        }
    }
    set.seed(NULL)

    print(c(cxy=mean(unlist(cxy)), cgy=mean(unlist(cgy)), cgx=mean(unlist(cgx))))

    rpt <- cbind(arg, do.call(rbind, rpt))
    pow(rpt)
}
