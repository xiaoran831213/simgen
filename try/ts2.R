tst2 <- function(N=5e2, M=2, L=4, a=0, b=1, d=0, e=1, times=1e3, seed=NULL, ...)
{
    arg <- get.arg(skp=c("seed", "times"))
    evt <- arg$evt %||% 1e-8
    psd <- arg$psd %||% evt / 10
    set.seed(seed)
    rhm <- ~G(L) + X(L) + U(L * 2) + D(M) | LD(G, a=8, b=.1) + CX(X, a=1, b=1) + CU(U, a=1, b=1) +
        ES(D, a=1, b=1) + .2 @ GX(G+X, a=1, b=1) + 0 @ GU(G+U, a=1, b=1)
    lhm <- Y(M) ~ a @ G + b @ X + d @ U:G
    rpt <- list()
    for(ii in seq(times))
    {
        flood(sim_dsg(rhm, N, psd=psd)$dat)
        G <- as.genotype(G)
        F <- fac(G)                    # genotype factorized
        G <- std(G, TRUE, FALSE)
        flood(sim_rsp(lhm)$dat)        # additive effect
        Q <- .p(G, 1:2, 1)             # genotype expanded
        Q <- Q[, asd(Q) > 0]           # remove degenerated terms
        Y <- std(Y + D * e)            # outcomes
        rG <- get_rsd(Y, X, G)         # residual, G as it is
        rF <- get_rsd(Y, X, F)         # residual, factorized G
        rQ <- get_rsd(Y, X, Q)         # residual, quadratic G
        LHS <- list(
            ## `LHS = Y^2` = .p(Y, 2:2, 1:1),
            ## `LHS = Y:Y` = .p(Y, 2:2, 2:2),
            `LHS = (Y-G-X)*` = .p(rG, 2:2, 1:2),
            `LHS = (Y-F-X)*` = .p(rF, 2:2, 1:2),
            `LHS = (Y-Q-X)*` = .p(rQ, 2:2, 1:2),
            `LHS = Y` = Y)
        LHS <- LHS[sapply(LHS, length) > 0]
        MDL <- list(
            `LHS ~ {G} + X` = lhs ~ {G} + X,
            `LHS ~ {F} + X` = lhs ~ {F} + X,
            `LHS ~ {Q} + X` = lhs ~ {Q} + X)
        CFG <- .d(lhs=names(LHS), mdl=names(MDL)[c(1:3, 1)])

        r <- list()
        for(j in seq(nrow(CFG)))
        {
            cfg <- CFG[j, ]
            lhs <- LHS[[cfg$lhs]] # left hand side
            mdl <- MDL[[cfg$mdl]] # model
            aso <- gwa_lm(mdl)    # association analysis
            lhc <- cor(aso$rsp)   # response correlation
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
