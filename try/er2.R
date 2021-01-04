#' Temporary Test
#'
#' @param N  # of samples
er2 <- function(N=5e2, M=2, a=0, times=1e3, ...)
{
    arg <- get.arg(skp=c("seed", "times"))
    evt <- arg$evt %||% 1e-8
    psd <- arg$psd %||% evt / 10
    set.seed(arg$seed)
    rpt <- list()
    for(ii in seq(times))
    {
        G <- kgp(N, 4, 1e-8, .99)   # get real genotype
        F <- fac(G)                     # factorized
        Y <- matrix(rnorm(N * M), N, M) # outcomes
        Y2 <- Y^2
        YY <- .p(Y, 2, 2)
        LHS <- list(
            `LHS = Y^2`      = Y2,
            `LHS = Y:Y`      = YY,
            `LHS = BCX(Y^2)` = BCX(Y2),
            `LHS = QRQ(Y:Y)` = ORQ(YY))
        G <- ORQ(G)
        F <- ORQ(F)
        MDL <- list(
            `LHS ~ {G}` = lhs ~ {G},
            `LHS ~ {F}` = lhs ~ {F})
        PLD <- lapply(MDL, gwa_pld)
        
        CFG <- .e(lhs=names(LHS), mdl=names(MDL))
        r <- list()
        for(j in seq(nrow(CFG)))
        {
            cfg <- CFG[j, ]
            lhs <- LHS[[cfg$lhs]] # left hand side
            mdl <- MDL[[cfg$mdl]] # model
            pld <- PLD[[cfg$mdl]] # partial LD
            aso <- gwa_lm(mdl)    # association analysis
            lhc <- cor(aso$rsp)   # response correlation
            . <- mdt(aso$zsc, pld, lhc, tol.egv=evt, ...)
            r <- r %c% .d(cfg, mtd="DOT", pvl=.$P, egv=.$L)
        }
        mcr <- mean(abs(cor(G)))
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
