tst6 <- function(N=5e2, M=2, a=0, b=1, d=0, e=1, times=1e3, ...)
{
    arg <- get.arg(skp=c("seed", "times"))
    evt <- arg$evt %||% 1e-8
    psd <- arg$psd %||% evt / 10
    set.seed(arg$seed)
    rhm <- ~G(5) + X(5) + U(9) | LD(G, a=6, b=1) + CX(X, a=1, b=1) + CU(U, a=1, b=1) +
        .5 @ GX(G+X, a=1, b=1) + .2 @ GU(G+U, a=1, b=1)
    lhm <- Y(M) ~ a @ G + b @ X + d @ U:G
    rpt <- list()
    for(ii in seq(times))
    {
        flood(sim_dsg(rhm, N, psd=psd)$dat)
        G <- fac(as.genotype(G))
        flood(sim_rsp(lhm)$dat)
        Y <- Y + matrix(rnorm(N * M), N, M) * e # outcomes
        Z <- lm(Y ~ X)$resid
        R <- lm(Y ~ X + G)$resid                # residual
        Y <- std(Y)
        LHS <- list(
            `LHS = Y`         = Y,
            `LHS = BCX(Z^2)`  = BCX(.p(Z, 2:2, 1:1)),
            `LHS = ORQ(Z:Z)`  = ORQ(.p(Z, 2:2, 2:2)),
            `LHS = Z^2`       = .p(Z, 2:2, 1:1),
            `LHS = Z:Z`       = .p(Z, 2:2, 2:2),
            `LHS = R^2 + R:R` = .p(R, 2:2, 1:2))
        MDL <- list(
            `LHS - X     ~ {G}` = lhs - X     ~ {G},
            `LHS - X - Y ~ {G}` = lhs - X - Y ~ {G})
        PLD <- lapply(MDL, gwa_pld)
        CFG <- rbind(
            .e(lhs=names(LHS), mdl=names(MDL)[1:1]))

        r <- list()
        for(j in seq(nrow(CFG)))
        {
            cfg <- CFG[j, ]
            lhs <- LHS[[cfg$lhs]] # left hand side
            mdl <- MDL[[cfg$mdl]] # model
            pld <- PLD[[cfg$mdl]] # partial LD
            aso <- gwa_lm(mdl)    # association analysis
            lhc <- cor(aso$rsp)   # response correlation
            . <- mdt(aso$zsc, pld, lhc, tol.egv=evt)
            r <- r %c% .d(cfg, mtd="DOT", pvl=.$P, egv=.$L)
            . <- mtq(aso$zsc, pld, lhc, tol.egv=evt)
            r <- r %c% .d(cfg, mtd="TQT", pvl=.$P, egv=.$L)
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
