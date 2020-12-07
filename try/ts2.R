#' Temporary Test
#'
#' @param N  # of samples
#' @param L  # of variants
#' @param P  # of casual variants
#' @param M  # of covariates
#' @param nsd SD of noise
tst2 <- function(N=2e2, es1=0, es2=0, eps=1, ydm=2, ydt=1,
                 cdt=5, evt=1e-4,
                 seed=NULL, times=3e3, ...)
{
    arg <- get.arg(skp=c("seed", "out", "times"))
    psd <- arg$psd %||% evt
    set.seed(seed)

    rpt <- list()
    for(ii in seq(times))
    {
        ## G: SNPs in high LD
        flood(sim_dsg(~G(5) | LD(G, a=3, b=1), psd=psd, N=N)$dat)
        G <- std(G)
        ## Y: correlated traits
        flood(sim_dsg(~Y(ydm) | CY(Y, a=1, b=1), psd=psd, N=N)$dat)
        Y <- std(Y)
        ## C: G and Y collide in covariates 
        C <- sim_rsp(C(cdt) ~ es1 @ G + es2 @ Y)
        C <- std(C) + matrix(rnorm(N * ncol(C)), N, ncol(C)) * eps

        ## residual of Y ~ C
        Z <- as.matrix(resid(lm(Y ~ C)))
        if(ydt==1)
            Z <- Z %*% nsp(cor(Z), NULL)$H
        if(ydt==2)
            Z <- pcs(Z)
        ## residual of G ~ C
        J <- as.matrix(resid(lm(G ~ C)))
        
        ## GWAS between nY and nX
        t0 <- gwa_lm(Y ~     {G}) # no covariate
        t1 <- gwa_lm(Y ~ C + {G}) # original t tests
        t2 <- gwa_lm(Z ~     {G}) # residual t tests
        t3 <- gwa_lm(Z ~     {J}) # double residual

        ## test statistics
        r <- list()
        r[['DT0']] <- mdt(t0$zsc, cor(G), cor(Y), tol.egv=evt) # controlled
        r[['DT1']] <- mdt(t1$zsc, cor(J), cor(Y), tol.egv=evt) # 
        r[['DT2']] <- mdt(t2$zsc, cov(J), cor(Z), tol.egv=evt) # 
        r[['DT3']] <- mdt(t3$zsc, cor(J), cor(Z), tol.egv=evt) # 

        p <- sapply(r, `[[`, 'P')
        l <- sapply(r, `[[`, 'L')
        rpt[[ii]] <- .d(itr=ii, mtd=names(r), pvl=p, egv=l)
        if(!(ii %% 10))
            cat(".", ii, ".", sep=""); if(!(ii %% 100)) cat("\r")
    }
    set.seed(NULL)

    rpt <- cbind(arg, do.call(rbind, rpt))
    rownames(rpt) <- NULL
    pow(rpt)
}
