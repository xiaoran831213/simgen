#' Temporary Test
#'
#' @param N  # of samples
#' @param M  # of traits
#' @param nsd SD of noise
tst0 <- function(N=2e2, M=2, a=0, e=1, evt=1e-4, seed=NULL, times=1e3, ...)
{
    arg <- get.arg(skp=c("seed", "out", "times"))
    psd <- arg$psd %||% evt
    set.seed(seed)

    rpt <- list()
    nG <- 5
    nU <- 0
    rhs <- ~G(nG + nU) | HL(G, a=1, b=1)
    for(ii in seq(times))
    {
        flood(sim_dsg(rhs, N, psd=psd)$dat)
        U <- G[, -seq(nG)] # unobserved
        G <- G[, +seq(nG)]  # SNPs
        flood(sim_rsp(Y(M) ~ a @ G)$dat)
        Y <- Y + matrix(rnorm(N * M, 0, e), N, M)
        ## GWAS
        Y <- std(Y)
        t0 <- gwa_lm(Y ~ {G}) # no covariate
        ## test statistics
        r <- list()
        r[['D0a']] <- mdt(t0$zsc, cor(G), cor(Y), tol.egv=evt) #
        r[['T0a']] <- mtq(t0$zsc, cor(G), cor(Y), tol.egv=evt) # controlled

        p <- sapply(r, `[[`, 'P')
        l <- sapply(r, `[[`, 'L')

        mcr <- mean(abs(cor(G, U)))
        rpt[[ii]] <- .d(itr=ii, mtd=names(r), pvl=p, egv=l, mcr=mcr)
        if(!(ii %% 10))
        {
            cat(".", ii, ".", sep=""); if(!(ii %% 100)) cat("\r")
        }
    }
    set.seed(NULL)

    rpt <- cbind(arg, do.call(rbind, rpt))
    pow(rpt)
}
