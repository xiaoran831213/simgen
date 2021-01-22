tst3 <- function(N=5e2, M=2, L=4, a=0, b=1, d=0, e=1, times=1e3, ...)
{
    arg <- get.arg(skp=c("seed", "times"))
    evt <- arg$evt %||% 1e-8
    psd <- arg$psd %||% evt / 10
    set.seed(arg$seed)
    lhm <- Y(M) ~ a @ G + b @ X + d @ U:G
    rhm <- ~U(10) + D(M) | CU(U, a=1, b=1) + ES(D, a=1, b=1)
    rpt <- list()
    for(ii in seq(times))
    {
        flood(sim_dsg(rhm, N, psd=psd)$dat)
        G <- kgp(N, L, kgp.con=1, kgp.psd=psd, kgp.ucr=.8)
        X <- std(kgp(N, 5, kgp.con=0, kgp.psd=psd, kgp.ucr=.8))
        F <- fac(G)                    # genotype factorized
        U <- std(get_rsd(U, G), TRUE, FALSE)
        flood(sim_rsp(lhm)$dat)        # additive effect
        Y <- std(Y + D * e)            # outcomes
        rG <- get_rsd(Y, X, G)         # residual, G as it is
        rF <- get_rsd(Y, X, F)         # residual, factorized G
        LHS <- list(
            ## `LHS = Y^2` = .p(Y, 2:2, 1:1),
            ## `LHS = Y:Y` = .p(Y, 2:2, 2:2),
            `LHS = (Y-G-X)*` = .p(rG, 2:2, 1:2),
            `LHS = (Y-F-X)*` = .p(rF, 2:2, 1:2),
            `LHS = Y` = Y)
        LHS <- LHS[sapply(LHS, length) > 0]
        MDL <- list(
            `LHS ~ {G} + X` = lhs ~ {G} + X,
            `LHS ~ {F} + X` = lhs ~ {F} + X)
        ## CFG <- .d(lhs=names(LHS), mdl=names(MDL)[c(1:3, 1)])
        CFG <- rbind(
            .d(lhs=names(LHS)[1:2], mdl=names(MDL)[1:2]),
            .e(lhs=names(LHS)[3:3], mdl=names(MDL)[1:1]))
        
        r <- list()
        for(j in seq(nrow(CFG)))
        {
            cfg <- CFG[j, ]
            lhs <- LHS[[cfg$lhs]] # left hand side
            mdl <- MDL[[cfg$mdl]] # model
            aso <- gwa_lm(mdl)    # association analysis
            lhc <- cor(aso$rsp)   # response correlation, adjusted
            pld <- cor(aso$gmx)   # partial LD, from partial genotype
            . <- mdt(aso$zsc, pld, lhc, tol.egv=evt)
            r <- r %c% .d(cfg, mtd="DOT", pvl=.$P, egv=.$L)
            . <- mtq(aso$zsc, pld, lhc, tol.egv=evt)
            r <- r %c% .d(cfg, mtd="TQT", pvl=.$P, egv=.$L)
        }
        rpt[[ii]] <- cbind(itr=ii, do.call(rbind, r))
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
