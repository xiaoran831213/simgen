#' minor allele frequencies
#'
#' @param g genotype matrix, N row samples and M column variants
#' @return MAF of M variants
maf <- function(g) {a <- colMeans(g, na.rm=TRUE) / 2; pmin(a, 1 - a)}

#' allele sandard deviation
#'
#' @param g genotype matrix, N row samples and M column variants
#' @return SD of M variants
asd <- function(g) apply(g, 2, sd)

#' get formula factors
#'
#' A convenience function equivalent to `attr(terms(f), "factors")`.
#' 
#' @param f the formula
#' @return factors in the formula
#' @noRd
ffc <- function(f) attr(terms(f), "factors")


#' Cached Real Genotype
#' 
#' A segment taken from from 1000 genome project.
#'
#' A continuous segment from chromosome 17  region q12, with no missed calls and
#' minor allele frequency no less than 0.05.
#'
#' @name c17q12
#' @noRd
if(!exists("c17"))
{
    c17 <- readRDS(system.file("extdata", '17q12.rds', package="simgen"))
    n17 <- nrow(c17)
    m17 <- ncol(c17)
    f17 <- maf(c17)
}


#' retain non-colinear variables
#'
#' Trim a correlation  matrix so the correlation among  remaining variables are
#' within specified thresholds.
#' @param ldm LD correlation matrix
#' @param lcr lower correlation threshold
#' @param ucr upper correlation threshold
#' @return index of variables to be retained.
ncv <- function(ldm, lcr=0, ucr=1, ...)
{
    r <- abs(ldm)
    r[upper.tri(r, TRUE)] <- lcr + (ucr - lcr) / 2
    b <- which(lcr <= r & r <= ucr, arr.ind=TRUE)
    apply(r, 2, function(.) all(lcr <= . & . <= ucr))
}

#' set random element to missing
#'
#' @param x true data
#' @param f frequency of missing
#' @param nan the NaN to use (def=NA).
#' @return the incomplete observation of \code{x}
set.nan <- function(x, f=.1, nan=NA)
{
    x[sample.int(length(x), length(x) * f)] <- nan
    x
}


#' Genotype from 1000 Genome Project
#'
#' Randomly draw a segment from the 1000 Genome Project.
#'
#' @param N number of samples
#' @param L number of SNPs
#' @param psd positive definite threshold
#' @param ucr upper correlation threshold
kgp <- function(N, L, psd=NULL, ucr=NULL, MAF=0.05, MSD=NULL, ...)
{
    psd <- psd %||% sqrt(.Machine$double.eps)
    ucr <- ucr %||% 0.99
    MAF <- MAF %||% 0.05
    MSD <- MSD %||% sqrt(2 * MAF * (1 - MAF)) * .8

    ## L variants on demand, reserve 2 * L
    P <- min(L * 5, m17)
    while(TRUE)
    {
        i <- sample.int(n17, N, N > n17)            # N
        j <- seq(sample(m17 - P, 1) + 1, l=P)       # P
        gmx <- c17[i, j]                            # N x P
        gmx <- gmx[, maf(gmx) >= MAF]               # minimum minor allele freq
        gmx <- gmx[, asd(gmx) >= MSD]               # minimum allele SD
        if(NCOL(gmx) < L)
        {
            P <- min(P + L - NCOL(gmx), m17)
            cat("P = ", P, "\n", sep="")
            next
        }
        ldm <- cor(gmx)                             # P x P
        
        ## drop  SNP  to  enforce  non-linearity;  if there  are  less  than  L
        ## remaining, try again with a bigger reserve.
        gmx <- gmx[, ncv(ldm, ucr=ucr, ...), drop=FALSE]
        if(NCOL(gmx) < L)
        {
            P <- min(P + L - NCOL(gmx), m17)
            cat("P = ", P, "\n", sep="")
            next
        }

        ## select L variants now
        j <- seq(sample(ncol(gmx) - L, 1) + 1, l=L)
        gmx <- gmx[, j, drop=FALSE]
        ldm <- cor(gmx)

        ## if min(eigenvalue) < (threshold) * max(eigenvalue), try again
        egv <- eigen(ldm, TRUE, TRUE)$values
        if (egv[L] < psd * egv[1L])
        {
            cat("Non-PSD!\n")
            next
        }
        break
    }
    gmx
}


#' Simulate variance covariance structure 
#'
#' @examples
#' rho <- sim_cor(8)
#' print(rho[1:5, 1:5])
#' summary(rho[lower.tri(rho)])
#' 
#' @noRd
sim_cor <- function(size, alpha=1, beta=alpha, tol.psd=NULL, ...)
{
    tol.psd <- tol.psd %||% sqrt(.Machine$double.eps)
    ## draw the lower triangle
    r <- matrix(0, size, size)
    r[lower.tri(r)] <- 2 * rbeta(size * (size - 1) / 2, alpha, beta) - 1
    ## make correlation
    r <- r + t(r)
    ## enforce PD and diagonal
    npd(r, 1, tol.psd=tol.psd)
}


## parameters for the left hand side (i.e., response variables)
pa0 <- function(size=1, ...)
{
    list(size=size, ...)
}

## parameters for the right hand, 1st moment
pa1 <- function(size=1, type=c("numeric", "genotype", "kgb"), ...)
{
    dict <- c(numeric=1, genotype=2, kgb=3)
    if(is.character(type))
        type <- match.arg(type)
    list(size=size, type=type, ...)
}

## parameters for the right hand, 2nd moment
pa2 <- function(vars, type=c(diag=1, rand=2, ones=3), alpha=NULL, beta=NULL, ...)
{
    dict <- eval(formals(sys.function())$type, parent.frame())
    errm <- paste(paste(names(dict), dict, sep="="), collapse=", ")
    v <- all.vars(substitute(vars))
    a <- alpha %||% 1
    b <- beta  %||% 1
    i <- Inf
    if(missing(type))
        type <- dict[1]
    else if(is.character(type))
        type <- match.arg(type, names(dict))
    else if(is.numeric(type) && type %in% seq_along(dict))
        type <- dict[type]
    else
        stop("must choose from (", errm, ")")
    pars <- switch(
        type,
        rand=list(alpha=a, beta=b),
        diag=list(alpha=i, beta=i),
        ones=list(alpha=1, beta=0))
    c(list(vars=v), pars)
}

#' Simulate design matrix
#'
#' @examples
#' mu <- "X(2, 'n') + G(3, 'g')"
#' sg <- "CX(X, 'r') + LD(G, a=8, b=1) + 0.1 @ XG(X+G, 'r')"
#' fm <- as.formula(paste("~", mu, "|", sg))
#' print(fm)
#'
#' res <- sim_dsg(fm, 1e3)
#'
#' ## data
#' print(round(head(res$raw), 3))  # raw data
#' print(round(head(res$dsg), 3))  # transformed
#' 
#' ## correlation / LD
#' print(res$vcs)       # suggested
#' sapply(res$dat, cor) # empirical
sim_dsg <- function(model, N=5e2, psd=NULL)
{
    ## env <- environment()
    env <- parent.frame()
    ## separate the mm1 and mm2
    . <- strsplit(as.character(model)[2], " *[|] *")[[1]]
    mm1 <- as.formula(paste("~", .[1])) # 1st moment - mean
    mm2 <- as.formula(paste("~", .[2])) # 2nd moment - covariance structure
    
    ## -------- treat 1st moment, mm1 --------
    ## extract coefficients, full terms, and names
    . <- rownames(attr(terms(mm1), "factors"))
    co1 <- as.numeric(sub("@.*$", "", sub("^[^@]*$", "1", .)))
    tm1 <- sub("^.*@", "", .)        # terms
    nm1 <- sub("[(].*[)]$", "", tm1) # names
    names(co1) <- nm1
    
    ## parse inner term parameters
    fn1 <- sapply(sub("^[^(]*", "pa1", tm1), str2lang) # parse
    names(fn1) <- nm1
    fn1 <- lapply(fn1, eval, env)                      #

    ## -------- treat the 2nd moment, mm2 --------
    ## extract coefficients
    . <- colnames(ffc(mm2))
    co2 <- as.numeric(sub("@.*$", "", sub("^[^@]*$", "1", .)))

    . <- rownames(ffc(mm2))
    rgx <- regexpr("[A-z][^(]*([(][^()]*[)])", .) # match
    tm2 <- regmatches(., rgx)                     # terms
    nm2 <- sub("[(].*[)]$", "", tm2)              # names
    names(co2) <- nm2

    ## parse inner term parameters
    fn2 <- sapply(sub("^[^(]*", "pa2", tm2), str2lang) # parse
    names(fn2) <- nm2
    fn2 <- lapply(fn2, eval, env)                      #

    ## -------- generate data --------
    vsz <- sapply(fn1, `[[`, 'size') # variable dimensions
    vcn <- rep(nm1, vsz)             # variable column names
    vcs <- list()                    # variance covariance struct
    ## raw matrix
    raw <- matrix(0, N, sum(vsz), dimnames=list(NULL, vcn))
    for(. in names(fn2))
    {
        msk <- vcn %in% fn2[[.]]$vars # column mask
        vcs[[.]] <- do.call(sim_cor, c(size=sum(msk), fn2[[.]], tol.psd=psd))
        colnames(vcs[[.]]) <- vcn[msk]
        rownames(vcs[[.]]) <- vcn[msk]
        raw[, msk] <- raw[, msk] + mvn(N, 0, co2[.] * vcs[[.]])
    }
    ## design matrix and data
    dsg <- raw
    dat <- list()
    for(. in names(fn1))
    {
        msk <- vcn == .
        fun <- paste0("as.", fn1[[.]]$type)
        dsg[, msk] <- do.call(fun, c(list(x=raw[, msk]), fn1[[.]])) * co1[.]
        dat[[.]] <- dsg[, msk, drop=FALSE]
    }

    ## return
    list(tm1=tm1, fn1=fn1, mm1=mm1, co1=co1,
         tm2=tm2, fn2=fn2, mm2=mm2, co2=co2,
         vcs=vcs, raw=raw, dsg=dsg, dat=dat)
}

#' Simulate response
#'
#' mu <- "A(2, 'g') + B(3, 'g') + C(2, 'n') + D(2, 'n')"
#' sg <- "Ph(A, 'd') + 2 @ LD(B, a=7, b=1) + .5 @ CV(B + C, a=1, b=6) + .2 @ JS(D, 'd')"
#' fm <- as.formula(paste("~", mu, "|", sg))
#' print(fm)
#' dt <- sim_dsg(fm, 1e3)
#' 
#' md <- Y(4) + Z(2) ~ .5 @ A + .p(B, 2) + .5 @ B : C + D | .5 @ A
#' md <- Y(4) + Z(2) ~ 0 @ A + .p(B, 2) + .5 @ B : C + D
#' md <- Y(4) + Z(2) ~ 0 | .5 @ A + 2 @ .p(B, 2) + 2 @ B : C + D
#' rs <- sim_rsp(md, dt$dat)
sim_rsp <- function(model, data=NULL, ...)
{
    if(is.null(data))
        env <- parent.frame()
    else
    {
        flood(data)
        env <- environment()
    }
    ## separate left and right hand size
    lhs <- frm_lhs(model)
    rhs <- frm_rhs(model)

    ## check variables for sample size
    . <- lapply(all.vars(rhs), get, env)
    N <- sapply(., nrow)
    N <- unlist(N[!sapply(N, is.null)])
    stopifnot(length(unique(N)) == 1)
    N <- unique(N)
    
    ## separate mm1 and mm2
    mm1 <- frm_lbs(rhs) # 1st moment - mean
    mm2 <- frm_rbs(rhs) # 2nd moment - variance
    
    ## -------- treat mm1 --------
    if(frm_is0(mm1))
    {
        dt1 <- matrix(0, N, 0)
        sz1 <- numeric()
        co1 <- numeric()
    }
    else if(frm_is1(mm1))
    {
        dt1 <- matrix(1, N, 1, dimnames=list(NULL, "It"))
        sz1 <- list(Int=1)
        co1 <- 1
    }
    else
    {
        ## extract coefficients
        . <- colnames(attr(terms(mm1), "factors"))
        co1 <- sub("@.*$", "", sub("^[^@]*$", "1", .))
        co1 <- lapply(lapply(co1, str2lang), eval, parent.frame())
        tm1 <- sub("^.*@", "", .) # terms
        ## get data
        . <- sapply(paste0("~", tm1), as.formula, env=env)
        dt1 <- lapply(., frm_mtx, simplify = FALSE)
        sz1 <- lapply(dt1, ncol)  # size(s)
        dt1 <- do.call(cbind, dt1)
    }
    
    ## -------- treat mm2 --------
    ## extract coefficients
    if(frm_is0(mm2))
    {
        dt2 <- matrix(0, N, 0)
        sz2 <- numeric()
        co2 <- numeric()
    }
    else if(frm_is1(mm2))
    {
        dt2 <- matrix(1, N, 1, dimnames=list(NULL, "It"))
        sz2 <- list(Int=1)
        co2 <- 1
    }
    else
    {
        . <- colnames(attr(terms(mm2), "factors"))
        co2 <- as.list(as.numeric(sub("@.*$", "", sub("^[^@]*$", "1", .))))
        tm2 <- sub("^.*@", "", .) # terms
        ## get data
        . <- sapply(paste("~", tm2), as.formula, env=env)
        dt2 <- sapply(., frm_mtx, simplify = FALSE)
        sz2 <- lapply(dt2, ncol) # size(s)
        dt2 <- do.call(cBind, dt2)
    }
    
    ## -------- treat response --------
    . <- rownames(attr(terms(lhs), "factors"))
    tm0 <- sub("^.*@", "", .)        # terms
    nm0 <- sub("[(].*[)]$", "", tm0) # names
    ## coefficients
    co0 <- as.numeric(sub("@.*$", "", sub("^[^@]*$", "1", .)))
    names(co0) <- nm0
    ## parse inner term parameters
    fn0 <- sapply(sub("^[^(]*", "pa0", tm0), str2lang) # parse
    names(fn0) <- nm0
    fn0 <- lapply(fn0, eval, env)                      #
    sz0 <- sapply(fn0, `[[`, 'size')                   # sizes
    cn0 <- rep(nm0, sz0)                               # column names
    ## generate
    .r0 <- function(n, s) rnorm(n, 0, s)
    ef1 <- lapply(cn0, function(.) unlist(mapply(.r0, sz1, co1, SIMPLIFY=FALSE)))
    ef2 <- lapply(cn0, function(.) unlist(mapply(.r0, sz2, co2, SIMPLIFY=FALSE)))
    rsp <- mapply(function(e1=NULL, e2=NULL)
    {
        m <- if(is.null(e1)) 0 else dt1 %*% e1
        v <- if(is.null(e2)) 1 else log(1 + exp(dt2 %*% e2))
        rnorm(N, m, sqrt(v))
    },
    ef1, ef2)
    colnames(rsp) <- cn0

    ## separated data
    dat <- sapply(nm0, function(.) rsp[, cn0 == .], simplify = FALSE)
    
    list(tm0=tm0, co0=co0, sz0=sz0, fn0=fn0,
         mm1=mm1, co1=co1, sz1=sz1, ef1=ef1,
         mm2=mm2, co2=co2, sz2=sz2, ef2=ef2,
         dat=dat, rsp=rsp)
}


#' Discretizing dosage by Hardy-Weinberg Equilibrium
#'
#' @param x matrix of allele dosages in continuous scale
#' @param m vector of allele frequencies
#' @return allele dosage in  {0, 1, 2}  with frequency {m^2, 2*m*(1-m), (1-m)^2}
#' @noRd
as.genotype <- function(x, m=NULL, ...)
{
    if (is.null(m))
    {
        m <- f17[sample(m17 - ncol(x), 1) + seq(ncol(x))]
    }
    for(j in seq(ncol(x)))
    {
        m2 <- m[j] * m[j]           # HWE freq of 0 allele
        mn <- 2 * m[j] * (1 - m[j]) # HWE freq of 1 allele
        qt <- quantile(x[, j], c(m2, m2 + mn))
        x[, j] <- 2 - (x[, j] > qt[1]) - (x[, j] > qt[2])
    }
    x
}
