#' Temporary Test
#'
#' @param N  # of samples
tst6 <- function(N=5e2, M=2, a=0, b=1, d=0, e=1, evt=1e-4, times=1e3, ...)
{
    arg <- get.arg(skp=c("seed", "times"))
    psd <- arg$psd %||% evt
    set.seed(arg$seed)
    rhm <- ~G(5) + X(5) + U(9) | LD(G, a=5, b=1) + CX(X, a=1, b=1) +
        CU(U, a=1, b=1)
    ## + 0.5 @ GX(G+X, a=1, b=1) + 0.1 @ GU(G+U, a=a, b=1)
    lhm <- Y(M) ~ a @ G + b @ X + d @ U:G
    rpt <- list()
    for(ii in seq(times))
    {
        flood(sim_dsg(rhm, N, psd=psd)$dat)
        X <- X + G %*% matrix(rnorm(ncol(G) * ncol(X)), ncol(G), ncol(X)) * .2
        U <- U + G %*% matrix(rnorm(ncol(G) * ncol(U)), ncol(G), ncol(U)) * .1
        flood(sim_rsp(lhm)$dat)
        
        Y <- Y + matrix(rnorm(N * M), N, M) * e
        LHS <- list(
            `LHS = Y`               = Y,
            ## `LHS = Y + Y^2`         = .p(Y, 1:2, 1),
            ## `LHS = Y + Y:Y`         = cbind(Y, .p(Y, 2, 2)),
            `LHS = Y + Y^2 + Y:Y`   = .p(Y, 1:2, 1:2),
            ## `LHS = Y^2`             = .p(Y, 2, 1),
            ## `LHS = Y:Y`             = .p(Y, 2, 2))
            `LHS = Y^2 + Y:Y`       = .p(Y, 2:2, 1:2))
        MDL <- list(
            ## `LHS ~ X   +   {G}` = lhs ~ X   +   {G},
            `LHS - X     ~ {G}` = lhs - X   ~   {G},
            ## `LSH ~ X + Y + {G}` = lhs ~ X + Y + {G},
            `LHS - X - Y ~ {G}` = lhs - X - Y ~ {G})
        PLD <- lapply(MDL, gwa_pld)

        .ex <- function(...) expand.grid(..., stringsAsFactors = FALSE)
        CFG <- rbind(.ex(lhs=names(LHS)[1:2], mdl=names(MDL)[1:1]),
                     .ex(lhs=names(LHS)[3:3], mdl=names(MDL)[2:2]))
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
