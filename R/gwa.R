gwa_frm <- function(model, ...)
{
    ## separate left and right side
    lh <- frm_lhs(model, 1)
    rh <- frm_rhs(model, 1)
    en <- environment(model)
    
    ## from the r-hand side, separate genotype term(s)
    pt <- "[{][^{}]*[}]"
    gt <- regmatches(rh, gregexpr(pt, rh), invert=0)[[1]]
    gt <- paste(collapse="+", gsub("^[{]|[}]$", "", gt[gt != ""]))

    ## remainig r-hand side, no genotype
    rh <- regmatches(rh, gregexpr(pt, rh), invert=1)[[1]]
    rh <- gsub("^[+]|[+]$", "", rh)
    rh <- paste(collapse="+", rh[rh != ""])
    
    ## from the l-hand side, separate residual regresion term(s)
    pt <- "-[^-+]*"
    rr <- regmatches(lh, gregexpr(pt, lh), invert=0)[[1]]
    lh <- regmatches(lh, gregexpr(pt, lh), invert=1)[[1]]
    rr <- paste(collapse="+", gsub("-", "", rr[rr != ""]))
    lh <- paste(collapse="+", gsub("^[+]|[+]$", "", lh[lh != ""]))
    
    ## return
    lh <- if(lh=="") NULL else as.formula(paste("~", lh), en)
    rh <- if(rh=="") NULL else as.formula(paste("~", rh), en)
    rr <- if(rr=="") NULL else as.formula(paste("~", rr), en)
    gt <- if(gt=="") NULL else as.formula(paste("~", gt), en)

    list(lh=lh, rh=rh, gt=gt, rr=rr)
}


#' GWAS linear model
#' @examples
#'
#' set.seed(1245)
#' N <- 100                          # sample size
#' D <- mvn(N, 0, cmx(9))            # 9 predictors in total
#' Z <- D[, 5:6]                     # 2 covariates observed
#' X <- hwe(D[, 1:4])                # 4 g-variants observed
#' D[, 1:4] <- X
#' U <- D[, 7:9]                     # 3 predictors unseen
#' D <- std(D)
#' 
#' b <- matrix(rnorm(9 * 3), 9, 3)   # effects
#' Y <- std(D %*% b)                 # 
#' Y <- Y + mvn(N, 0, cmx(3))        # noise
#'
#' Y <- std(Y)
#' X <- std(X)
#' Z <- std(Z)
#'
#' y1 <- gwa_rsp(Y - Z ~ {X})        # residual response
#' y2 <- gwa_rsp(Y ~ {X} + Z)        # original response
#'
#' round(cor(y1, Z), 2)              # expecting zeros
#' round(cor(y2, Z), 2)              # expecting non-zeros
#'
#' round(colMeans(y1), 2)            # expecting zeros
#' round(colMeans(y2), 2)            # expecting zeros
gwa_lm <- function(model, ...)
{
    ## extract the rsponse, genotype, and covariate
    flood(gwa_frm(model, ...))

    rsp <- frm_mtx(lh)      # response
    gmx <- frm_mtx(gt, ...) # genotype
    cvr <- frm_mtx(rh, ...) # covariate
    ror <- frm_mtx(rr)      # residual regression
    if(is.null(rr))
        rsp <- resid(lm(rsp ~   1))
    else
        rsp <- resid(lm(rsp ~ ror))

    ## model
    mdl <- if(is.null(cvr)) y ~ g else y ~ g + cvr

    ## z-scores
    zsc <- matrix(0, ncol(gmx), ncol(rsp))
    for(i in seq(ncol(gmx)))
    {
        for(j in seq(ncol(rsp)))
        {
            g <- gmx[, i]
            y <- rsp[, j]
            m <- lm(mdl)
            zsc[i, j] <- summary(m)$coef[2, 3]
        }
    }
    zsc
}


gwa_cst <- function(md, ...)
{
    ## extract the rsponse, genotype, and covariate
    flood(gwa_frm(md, ...))

    ## get correlation among response variables
    rsp <- frm_mtx(lh)      # response
    ror <- frm_mtx(rr)      # residual regression
    if(is.null(ror))
        rsp <- resid(lm(rsp ~   1))
    else
        rsp <- resid(lm(rsp ~ ror))
    rsp <- cor(rsp)

    ## get correlation among genotype, controlled for covariates
    cvr <- frm_mtx(rh, ...) # covariate
    gmx <- frm_mtx(gt, ...) # genotype
    if(is.null(cvr))
        gmx <- cor(gmx)
    else
        gmx <- ginv(ginv(cor(cbind(gmx, cvr)))[1:ncol(gmx), 1:ncol(gmx)])

    list(gmx=gmx, rsp=rsp)
}
