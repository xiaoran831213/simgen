#' Compute Polynomials
#'
#' @param x matrix of basic terms.
#' @param d degrees to compute.
#' @param t number of terms allowed.
#' @param o use orthoganal transformation (def=FALSE).
#' @return polynomial terms expanded from the given basic terms.
#'
#' @export
.p <- function(x, d=NULL, t=NULL, o=FALSE)
{
    if(!is.matrix(x))
        dim(x) <- c(length(x), 1L)
    d <- d %||% 1
    t <- t %||% seq.int(max(d))

    ## get polynomial
    x <- stats::poly(x, degree=max(d), raw=!o, simple=TRUE)

    ## select by degree
    x <- x[, attr(x, 'degree') %in% d, drop=FALSE]

    ## select by terms
    . <- rowSums(do.call(rbind, strsplit(colnames(x), "[.]")) != "0")
    x <- x[, . %in% t, drop=FALSE]

    ## clean up
    attr(x, 'class') <- NULL
    attr(x, 'degree') <- NULL

    x
}

#' set random element to missing
#'
#' @param x true data
#' @param f frequency of missing
#' @param nan the NaN to use (def=NA).
#' @return the incomplete observation of \code{x}
#' @noRd
set.nan <- function(x, f=.1, nan=NA)
{
    x[sample.int(length(x), length(x) * f)] <- nan
    x
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
    n <- size * (size - 1) / 2
    r[lower.tri(r)] <- 2 * rbeta(n, alpha, beta) - 1
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
#' @param model formula to specify the design matrices
#' @param N sample size
#' @param psd threshold of positive definite correlation
#' @return list of parsed model and simulated data
#' 
#' @examples
#' ## print(sim_dsg)
#' ## d <- sim_dsg(~X(2) + G(3, 'g') | CX(X, 'r') + LD(G, a=8, b=1) + 0.1 @ XG(X+G, 'r'), 1e3)
#'
#' ## print(round(head(d$raw), 3))  # raw data
#' ## print(round(head(d$dsg), 3))  # transformed
#' 
#' ## print(d$vcs)       # suggested
#' ## sapply(d$dat, cor) # empirical
sim_dsg <- function(model, N, psd=NULL)
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
#' @param model the data generation model in R formula.
#' @param data a list or environment with data objects.
#' @param ... additional arguments (absorb and ignore)
#' @return a list of parsed model and response variable.
#'
#' @examples
#' ## generate design matrix
#' ## d <- sim_dsg(~A(3, 'g') + B(2, 'n') | LD(A, a=8, b=1) + CV(B, 'r') + .1 @ CF(B+C, 'r'), 1e3)
#'
#' ## generate 5 responses
#' ## r <- sim_rsp(Y(3) + Z(2) ~ .5 @ A + .p(B, 2) + .5 @ A:B | .5 @ A, dt$dat)
#' ## genotype A affects both mean and variance (i.e., vQTL)
#'
#' ## preview
#' ## round(head(r$rsp), 3)
#' @export
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
        dt2 <- do.call(cbind, dt2)
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
