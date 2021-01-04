tst4 <- function(N=5e2, M=2, L=4, a=0, b=0, d=0, e=1, times=1e3, ...)
{
    arg <- get.arg(skp=c("seed", "times"))
    evt <- arg$evt %||% 1e-8
    psd <- arg$psd %||% evt / 10
    set.seed(arg$seed)
    rhm <- ~G(L) + X(L) + U(L) | LD(G, a=1, b=1) + CX(X, a=1, b=1) + CU(U, a=1, b=1)
    lhm <- Y(M) ~ a @ G + b @ X + d @ U:G
    rpt <- list()
    for(ii in seq(times))
    {
        flood(sim_dsg(rhm, N, psd=psd)$dat)
        mcr <- mean(abs(cor(G, X))) + mean(abs(cor(G, U)))
        G <- as.genotype(G)


        flood(sim_rsp(lhm)$dat)
        Y <- Y + matrix(rnorm(N * M), N, M) * e # outcomes
        G <- ORQ(fac(G))                        # as factor, normalized
        R <- get_rsd(Y, X, G)                   # residual
        Y <- std(Y)                             # centered
        LHS <- list(
            `LHS = BCX(Y^2)`  = BCX(.p(Y, 2:2, 1:1)),
            `LHS = ORQ(Y:Y)`  = ORQ(.p(Y, 2:2, 2:2)),
            `LHS = Y:Y`       = .p(Y, 2:2, 2:2),
            `LHS = Y^2`       = .p(Y, 2:2, 1:1),
            `LHS = R^2 + R:R` = .p(R, 2:2, 1:2))
        MDL <- list(`LHS ~ {G} + X` = lhs ~ {G} + X)
        
        CFG <- .e(lhs=names(LHS), mdl=names(MDL))
        r <- list()
        for(j in seq(nrow(CFG)))
        {
            cfg <- CFG[j, ]
            lhs <- LHS[[cfg$lhs]] # left hand side
            mdl <- MDL[[cfg$mdl]] # model
            aso <- gwa_lm(mdl)    # association analysis
            lhc <- cor(aso$rsp)   # response correlation
            pld <- cor(aso$gmx)
            . <- mdt(aso$zsc, pld, lhc, tol.egv=evt)
            r <- r %c% .d(cfg, mtd="DOT", pvl=.$P, egv=.$L)
            . <- mtq(aso$zsc, pld, lhc, tol.egv=evt)
            r <- r %c% .d(cfg, mtd="TQT", pvl=.$P, egv=.$L)
        }
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
