tst0 <- function(N=5e2, M=2, L=4, a=0, b=1, d=0, e=1, times=1e3, ...)
{
    arg <- get.arg(skp=c("seed", "times"))
    evt <- arg$evt %||% 1e-8
    psd <- arg$psd %||% evt / 10
    set.seed(arg$seed)
    rhm <- ~G(L) + X(L) + U(2 * L) | LD(G, a=1, b=1) + CX(X, a=1, b=1) + CU(U, a=1, b=1)
    lhm <- Y(M) ~ a @ G + b @ X + d @ U:G
    rpt <- list()
    for(ii in seq(times))
    {
        flood(sim_dsg(rhm, N, psd=psd)$dat)
        flood(sim_rsp(lhm)$dat)
        Y <- Y + matrix(rnorm(N * M), N, M) * e # outcomes
        Y <- std(Y)
        LHS <- list(`LHS = Y`             = Y)
        MDL <- list(`LHS ~ {G} + X` = lhs ~ {G} + X)
        PLD <- lapply(MDL, gwa_pld)
        CFG <- .e(lhs=names(LHS), mdl=names(MDL))
        r <- list()
        for(j in seq(nrow(CFG)))
        {
            cfg <- CFG[j, ]
            lhs <- LHS[[cfg$lhs]] # left hand side
            mdl <- MDL[[cfg$mdl]] # model
            ## pld <- PLD[[cfg$mdl]] # partial LD
            aso <- gwa_lm(mdl)    # association analysis
            lhc <- cor(aso$rsp)   # response correlation, adjusted
            pld <- cor(aso$gmx)
            . <- mdt(aso$zsc, pld, lhc, tol.egv=evt)
            r <- r %c% .d(cfg, mtd="DOT", pvl=.$P, egv=.$L)
            ## . <- mtq(aso$zsc, pld, lhc, tol.egv=evt)
            ## r <- r %c% .d(cfg, mtd="TQT", pvl=.$P, egv=.$L)
        }
        mcr <- mean(abs(cor(G, X))) + mean(abs(cor(G, U)))
        rpt[[ii]] <- cbind(itr=ii, do.call(rbind, r), mcr=mcr)
        if(!(ii %% 10))
            cat(".", ii, ".", sep=""); if(!(ii %% 100)) cat("\r")
    }
    set.seed(NULL)
    rpt <- cbind(arg, do.call(rbind, rpt))
    
    rpt <- within(rpt,
    {
        lhs <- factor(lhs, names(LHS))      # residuals
        mdl <- factor(mdl, names(MDL))      # left hand side
        mtd <- factor(mtd, c("DOT", "TQT")) # which test
    })
    pow(rpt)
}
