frm_fct <- function(f) attr(terms(f), "factors")

frm_lhs <- function(f, ret=0)
{
    e <- environment(f)
    f <- as.character(f)[-1]
    if(length(f) < 2)
        return(NULL)

    f <- gsub("[\n ]*", "", f[1])
    if(ret == 0)
        f <- as.formula(paste("~", f), e)
    f
}

frm_rhs <- function(f, ret=0)
{
    e <- environment(f)
    f <- as.character(f)[-1]
    f <- gsub("[\n ]*", "", f[length(f)]) # right
    if(ret == 0)
        f <- as.formula(paste("~", f), e)
    f
}

frm_lbs <- function(f, ret=0)
{
    e <- environment(f)
    f <- as.character(f)
    f <- f[length(f)]
    f <- strsplit(f, "[|]")[[1]][1]
    if(ret == 0)
        f <- as.formula(paste("~", f), e)
    f
}

frm_rbs <- function(f, ret=0)
{
    e <- environment(f)
    f <- as.character(f)
    f <- f[length(f)]
    f <- strsplit(f, "[|]")[[1]]
    if(length(f) < 2)
        f <- 0
    else
        f <- f[2]
    if(ret == 0)
        f <- as.formula(paste("~", f), e)
    f
}

frm_bar <- function(rhs, ...)
{
    e <- environment(rhs)
    f <- frm_rhs(rhs, 1)
    f <- strsplit(f, "[|]")[[1]]
    f <- paste("~", f)
    f <- sapply(f, as.formula, env = e)
    unname(f)
}

frm_is0 <- function(f, ...)
{
    attr(terms(f), "intercept") == 0 && length(all.vars(f)) == 0
}

frm_is1 <- function(f, ...)
{
    attr(terms(f), "intercept") == 1 && length(all.vars(f)) == 0
}

frm_cat <-function(..l, ..r, ..o="~")
{
    ..f <- paste(as.character(..l)[-1], ..o, as.character(..r)[-1])
    if(..o != "~")
        ..f <- paste("~", ..f)

    flood(environment(..l))
    flood(environment(..r))

    as.formula(..f)
}

frm_mtx <- function(f, ...)
{
    if(is.null(f))
        return(NULL)
    x <- model.matrix(update.formula(f, ~ . - 1), ...)
    attr(x, "assign") <- NULL
    x
}
