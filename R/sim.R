#' Simulate variance covariance structure 
#'
#' @examples
#' rho <- cmx(500)
#' print(rho[1:5, 1:5])
#' summary(rho[lower.tri(rho)])
#' 
#' @noRd
sim_vcs <- function(L, alpha=0.9, beta=alpha, ...)
{
    ## correlation
    R <- matrix(0, L, L)

    ## draw the lower triangle
    R[lower.tri(R)] <- 2 * rbeta(L * (L - 1) / 2, alpha, beta) - 1
    R <- R + t(R)
    diag(R) <- 1

    ## boost
    u <- tcrossprod(2 * rbeta(L, log(L) * 2, log(L) * 2) - 1) * L
    R <- cov2cor(R + u)

    ## force PD
    R <- fpd(R)
    cov2cor(R)
}

sim_cor <- function(size, alpha=1, beta=alpha, center=0, debug=0, ...)
{
    ## correlation
    r <- matrix(0, size, size)
    ## draw the lower triangle
    if(center)
        r[lower.tri(r)] <- 2 * rbeta(size * (size - 1) / 2, alpha, beta) - 1
    else
        r[lower.tri(r)] <- rbeta(size * (size - 1) / 2, alpha, beta)
    ## make correlation
    r <- r + t(r)
    diag(r) <- 1
    r <- cov2cor(fpd(r))
    
    if(debug)
    {
        with(eigen(r, TRUE),
        {
            cat("eigen values =\n")
            print(values)
            cat("square roots =\n")
            print(vectors %*% (sqrt(values) * t(vectors)))
            cat("summary of lower triangle:\n")
            print(summary(r[lower.tri(r)]))
        })
    }
    r
}

## make dosage genotype
#' Simulation formula
#'
#' @examples
#' f <- Y(3, d=2) | Z(2) * U ~ A(3, "g") + .p(B(5, "g"), 2, 1) + v(.9) | E(5) * F(2) + v(.1)
#' r <- sim_frm(f)
sim_frm <- function(model, ...)
{
    env <- environment(model)
    ## separate left and right side, divided by bars "|".
    lhs <- frm_bar(frm_lhs(model))
    rhs <- frm_bar(frm_rhs(model))

    ## get terms
    lht <- lapply(lapply(sapply(lhs, terms), attr, "factors"), rownames)
    rht <- lapply(lapply(sapply(rhs, terms), attr, "factors"), rownames)
    
    ## variable names
    ## lhv <- lapply(lht, sub, pattern="[(].*[)]", replacement="")
    ## rhv <- lapply(rht, sub, pattern="[(].*[)]", replacement="")
    
    ## ## data generation parameters
    ## tmp <- lapply(lht, gsub, pattern="^[^(]*", replacement="")
    ## lhp <- lapply(tmp, gsub, pattern="^$", replacement="()")
    ## tmp <- lapply(rht, gsub, pattern="^[^(]*", replacement="")
    ## rhp <- lapply(tmp, gsub, pattern="^$", replacement="()")
    
    ## ## data generation calls
    ## lhc <- mapply(lhv, lhp, FUN=function(v, p) paste0(v, " <- sim_par", p))
    ## rhc <- mapply(rhv, rhp, FUN=function(v, p) paste0(v, " <- sim_par", p))
    
    ## lhc <- lapply(lhc, function(blk)
    ## {
    ##     blk <- lapply(blk, function(exe)
    ##     {
    ##         eval(str2lang(exe))
    ##     })
    ##     blk
    ## })
    ## lhc
    rht
}

## parameters for the 1st moment
pa1 <- function(size=1, type=c("numeric", "integer", "genotype"), ...)
{
    dict <- c(numeric=1, genotype=2)
    if(is.character(type))
        type <- match.arg(type)
    list(size=size, type=type, ...)
}

pa2 <- function(vars, type=c("rand", "diag", "ones"), alpha=NULL, beta=NULL, ...)
{
    v <- all.vars(substitute(vars))
    a <- alpha %||% 1
    b <- beta  %||% 1
    i <- Inf
    if(is.character(type))
        type <- match.arg(type)

    pars <- switch(
        type,
        rand=list(alpha=a, beta=b),
        diag=list(alpha=0, beta=1),
        ones=list(alpha=1, beta=0))
    c(list(vars=v), pars)
}

ex1 <- function(N=20)
{
    mu <- "A(2, 'g') + B(3, 'g') + C(2, 'n') + D(2, 'i')"
    ## sg <- "Phi(A, 'd') + 2 @ LD(B, 'r', a=.2, b=.2) + .5 @ CV(B + C, 'r', a=4, b=4)"
    sg <- "Phi(A, 'd') + 2 @ LD(B, 'r', a=.2, b=.2) + .5 @ CV(C, 'r', a=4, b=4)"
    fm <- paste(mu, "|", sg)
    frm <- as.formula(paste("~", fm))
    env <- environment()

    ## separate the mm1 and mm2
    . <- strsplit(fm, " *[|] *")[[1]]
    mm1 <- as.formula(paste("~", .[1])) # 1st moment - mean
    mm2 <- as.formula(paste("~", .[2])) # 2nd moment - covariance structure
    
    ## -------- treat 1st moment, mm1 --------
    ## extract coefficients
    . <- colnames(attr(terms(mm1), "factors"))
    co1 <- as.numeric(sub("@.*$", "", sub("^[^@]*$", "1", .)))
    mm1 <- as.formula(paste("~", paste(sub("^.*@", "", .), collapse = "+")))
    ## extract inner terms
    . <- rownames(attr(terms(mm1), "factors"))
    rgx <- regexpr("[A-z][^(]*([(][^()]*[)])", .) # match
    tm1 <- regmatches(., rgx)                     # terms
    nm1 <- sub("[(].*[)]$", "", tm1)              # names

    ## substitute inner terms with names
    . <- as.character(mm1)[2]
    rgx <- gregexpr("[A-z][^(]*([(][^()]*[)])", .) # match
    fc1 <- regmatches(., rgx)                      # terms
    regmatches(., rgx) <- lapply(fc1, sub, pattern="[(].*[)]$", replacement="")
    mm1 <- as.formula(paste("~", .))
    
    ## parse inner term parameters
    fn1 <- sapply(sub("^[^(]*", "pa1", tm1), str2lang) # parse
    names(fn1) <- nm1
    dt1 <- lapply(fn1, eval, env)                      #
    
    ## -------- treat the 2nd moment, mm2 --------
    ## extract coefficients
    . <- colnames(attr(terms(mm2), "factors"))
    co2 <- as.numeric(sub("@.*$", "", sub("^[^@]*$", "1", .)))
    mm2 <- as.formula(paste("~", paste(sub("^.*@", "", .), collapse = "+")))

    . <- rownames(attr(terms(mm2), "factors"))
    rgx <- regexpr("[A-z][^(]*([(][^()]*[)])", .) # match
    tm2 <- regmatches(., rgx)                     # terms
    nm2 <- sub("[(].*[)]$", "", tm2)              # names
    
    ## substitute inner terms with names
    . <- as.character(mm2)[2]
    rgx <- gregexpr("[A-z][^(]*([(][^()]*[)])", .) # match
    fc2 <- regmatches(., rgx)                      # terms
    regmatches(., rgx) <- lapply(fc2, sub, pa="[(].*[)]$", re="")
    mm2 <- as.formula(paste("~", .))

    ## parse inner term parameters
    fn2 <- sapply(sub("^[^(]*", "pa2", tm2), str2lang) # parse
    names(fn2) <- nm2
    dt2 <- lapply(fn2, eval, env)                      #

    ## -------- generate data --------
    dm1 <- sapply(dt1, `[[`, 'size') # dimensions
    mk1 <- rep(nm1, dm1)             # masks
    
    ## variance and covariance structure (vcs)
    vcs <- lapply(dt2, function(par)
    {
        msk <- ! mk1 %in% par$vars
        par <- c(size=sum(dm1), par)
        ret <- do.call(sim_cor, par)
        ret[msk, msk] <- 0
        ret
    })
    sgm <- Reduce("+", mapply("*", co2, vcs, SIMPLIFY = FALSE))
    ## raw data
    raw <- mvn(N, 0, sgm)
    
    list(tm1=tm1, fn1=fn1, dt1=dt1, mm1=mm1, co1=co1,
         tm2=tm2, fn2=fn2, dt2=dt2, mm2=mm2, co2=co2,
         vcs=vcs, sgm=sgm, raw=raw)
}



#' Discretizing dosage by Hardy-Weinberg Equilibrium
#'
#' @param x matrix of allele dosages in continuous scale
#' @param m vector of allele frequencies
#' @return allele dosage in  {0, 1, 2}  with frequency {m^2, 2*m*(1-m), (1-m)^2}
#' @noRd
hwe <- function(x, m=NULL, ...)
{
    if (is.null(m))
        m <- runif(ncol(x), 0.05, 0.45)
    for(j in seq(ncol(x)))
    {
        m2 <- m[j] * m[j]           # HWE freq of 0 allele
        mn <- 2 * m[j] * (1 - m[j]) # HWE freq of 1 allele
        qt <- quantile(x[, j], c(m2, m2 + mn))
        x[, j] <- 2 - (x[, j] > qt[1]) - (x[, j] > qt[2])
    }
    x
}

#' Make the right side
#'
#' @param N sample size
#' @param L number of g-variants
#' @param M number of covariates
#' @param P number of g-variants that is effective
#' @param e size of white noise
#' @return
#' a list,
#' \itemize{
#' }
sim <- function(nrep, expr, N=1e3, L=3, M=2, P=1, nsd=1, ...)
{
    cfg <- get.arg(skp=c('nrep', 'expr', 'maf'))
    bug <- cfg$bug %||% 0

    ## default expression: summarize the outcome Y
    if(missing(expr))
    {
        expr <- quote(
        {
            val <- summary(Y)
            key <- names(val)
            val <- round(as.vector(val), 3)
            .d(mtd="sm5", key=key, val=val)
        })
    }

    ## correlation of all variables
    C <- cmx(L + M, ...)

    res <- list()
    for(i in seq(nrep))
    {
        ## L variants and M covariates
        X <- mvn(N, 0, C)          

        ## variants discretized
        if(L > 0)
            X[, 1:L] <- hwe(X[, 1:L], ...)

        G <- X[, 0 + seq(1, l=L), drop=FALSE]
        X <- X[, L + seq(1, l=M), drop=FALSE]
        
        ## Y = noise + genetics + environment
        Y <- rnorm(N, 0, nsd) # noise
        if(L > 0 && P > 0)         # genetics
            Y <- Y + std(G[, sample(L, P), drop=FALSE] %*% rnorm(P))
        if(M > 0) # enviornment
            Y <- Y + std(X %*% rnorm(M))
        Y <- drop(std(Y))

        val <- eval(substitute(expr))
        res[[i]] <- .d(itr=i, val)
    }
    res <- cbind(cfg, do.call(rbind, res))
    
    ## return
    res
}
