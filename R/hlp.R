#' Multivariant Normal Samples
#'
#' A simplified copy of MASS::mvrnorm.
#'
#' @param N number of samples
#' @param M vector of mean,  rotated  to fit the  dimension of \code{V}
#' @param V matrix of covariance
#' @param drop TRUE to drop matrix of single sample to a vector
#' @return N x M matrix of samples
#' @noRD
mvn <- function (N=1, M=0, V=NULL, drop=TRUE)
{
    ## default V is for demonstration
    if(is.null(V))
        V = rbind(c(1, .5, .25), c(.5, 1, .5), c(.25, .5, 1))

    ## dimensionality
    D <- nrow(V)

    ## mean vector
    M <- drop(rep(M, length=D))

    ## eigen decomposition
    e <- eigen(V, symmetric = TRUE)
    d <- e$values
    U <- e$vectors
    s <- sqrt(pmax(d, 0))          # square root of V

    ## random values
    X <- matrix(rnorm(D * N), D, N)
    y <- t(M + U %*% (s * X))

    if(drop)
        y <- drop(y)
    y
}


#' Force Positive Definite
#' @noRD
fpd <- function(x, eps=NULL)
{
    if (is.null(eps))
        eps <- sqrt(.Machine$double.eps)

    e <- eigen(x, TRUE)
    d <- e$values
    u <- e$vectors

    d[d < eps] <- eps
    x <- u %*% (t(u) * d)
    x <- 0.5 * (x + t(x))
    x
}

min.evl <- function(x)
{
    min(eigen(x, TRUE, TRUE)$values)
}


#' Short for scale()
std <- function(x, c=TRUE, s=TRUE, ...)
{
    r <- scale(x, center=c, scale=s)
    attr(r, 'scaled:center') <- NULL
    attr(r, 'scaled:scale') <- NULL
    r
}

#' Compute Polynomials
#'
#' @param x matrix of basic terms
#' @param d degrees to compute
#' @param t number of terms allowed
#' @return polynomial terms expanded from the given basic terms.
.p <- function(x, d=NULL, t=NULL)
{
    if(!is.matrix(x))
        dim(x) <- c(length(x), 1L)
    d <- d %||% 1
    t <- t %||% seq.int(max(d))

    ## get polynomial
    x <- poly(x, degree=max(d), raw=TRUE, simple=TRUE)

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
