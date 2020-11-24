#' Temporary Test
#'
#' @param N  # of samples
#' @param L  # of variants
#' @param P  # of casual variants
#' @param M  # of covariates
#' @param nsd SD of noise
test <- function(N=2e2, es1=0, es2=0, eps=1, evt=1e-4, seed=NULL, times=3e3, ...)
{
    arg <- get.arg(skp=c("seed", "out", "times"))
    set.seed(seed)
    
    rpt <- list()
    for(ii in seq(times))
    {
        ## mdl <- ~F(1) + G(4) + C(5) | LD(G, a=6, b=1) + CC(C, a=1, b=1)
        mdl <- ~F(1) + G(4) + C(5) + U(10) | LD(F + G, a=2, b=1) + CC(C, a=.2, b=1) + UC(U, a=0, b=1)
        rhs <- sim_dsg(mdl, N)
        flood(rhs$dat)
        G <- std(G)
        C <- std(C)
        Y <- sim_rsp(Y(2) ~ es1 @ G + es2 @ G:U + C)
        Y <- Y + matrix(rnorm(length(Y)), nrow(Y), ncol(Y)) * eps
        Y <- std(Y)
        Y <- Y %*% nsp(cor(Y), NULL)$H
        
        ## cor among Y and X, the latter is controlled of Conf
        ## Cr <- gwa_cst(Y ~ C + {G})$gmx
        Cr <- scp(cor(cbind(G, C)), colnames(C))

        ## GWAS between nY and nX
        tt <- gwa_lm(Y ~ C + {G}) # original t tests
        rt <- gwa_lm(Y - C ~ {G}) # residual t tests

        ## test statistics
        r <- list()
        Cy <- diag(ncol(Y))
        r[['TQ1']] <- mtq(tt$zsc, cov2cor(Cr), Cy, tol.egv=evt)
        r[['DT1']] <- mdt(tt$zsc, cov2cor(Cr), Cy, tol.egv=evt)
        r[['TQ2']] <- mtq(rt$zsc, Cr, Cy, tol.egv=evt)
        r[['DT2']] <- mdt(rt$zsc, Cr, Cy, tol.egv=evt)

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

pow <- function(rpt)
{
    rpt <- subset(rpt, se=-itr)
    grp <- subset(rpt, se=-c(pvl, egv))
    rpt <- by(rpt, grp, function(g)
    {
        cfg <- subset(g, se=-c(pvl, egv))[1, ]
        pow <- with(g, mean(pvl <= 0.05))
        egv <- with(g, mean(egv))
        cbind(cfg, pow=pow, egv=egv, rep=nrow(g))
    })
    rpt <- do.call(rbind, rpt)
    rpt
}
