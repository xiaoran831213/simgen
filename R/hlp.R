#' Multivariant Normal Samples
#'
#' A simplified copy of MASS::mvrnorm.
#'
#' @param N number of samples
#' @param M vector of mean,  rotated  to fit the  dimension of \code{V}
#' @param V matrix of covariance
#' @param drop TRUE to drop matrix of single sample to a vector
#' @return N x M matrix of samples
#' @noRd
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


#' Near Positive Definite
#'
#' A simplified copy of Matrix::nearPD that project a squre matrix to near PD.
#'
#' @param x square matrix input
#' @param d the desired diagnal, using 1 forces the output to be a correlation matrix.
#' @param tol.egv eigen trimming threshold
#' @param tol.psd minimum eigenvalue
#' @param tol.cnv tolerence of convergence
#' @noRd
npd <- function (x, dg=NULL, tol.egv=1e-06, tol.cnv=1e-07, tol.psd=1e-08, max.itr=100)
{
    n <- ncol(x)
    if(!is.null(dg))
        dg <- rep(dg, length.out=n)

    ## enforce PD, trim components
    X_diff <- x * 0
    X_this <- x
    for(itr in seq(max.itr))
    {
        X_prev <- X_this
        X_adjt <- X_prev - X_diff

        e <- eigen(X_adjt, symmetric = TRUE)
        Q <- e$vectors
        D <- e$values
        . <- D > tol.egv * D[1]
        Q <- Q[, ., drop = FALSE]
        D <- D[.]
        X_this <- Q %*% (D * t(Q))
        ## changes due to trimming
        X_diff <- X_this - X_adjt

        ## enforce diagonal
        if(!is.null(dg))
            diag(X_this) <- dg

        ## total relative change
        if(norm(X_prev - X_this, 'I') / norm(X_prev, 'I') < tol.cnv)
            break
    }
    if (itr == max.itr)
        warning(gettextf("'nearPD()' did not converge in %d iterations", itr), domain=NA)

    ## enforece minimum eigen values
    e <- eigen(X_this, symmetric = TRUE)
    Q <- e$vectors
    D <- e$values
    eps <- tol.psd * abs(D[1])
    if (D[n] < eps)
    {
        D[D < eps] <- eps
        D0 <- diag(X_this)
        X_this <- Q %*% (D * t(Q))
        S <- sqrt(pmax(eps, D0) / diag(X_this))
        X_this[] <- outer(S, S) * X_this
    }
    if(!is.null(dg))
        diag(X_this) <- dg
    
    X_this
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
#' @param x matrix of basic terms.
#' @param d degrees to compute.
#' @param t number of terms allowed.
#' @param o use orthoganal transformation (def=FALSE).
#' @return polynomial terms expanded from the given basic terms.
.p <- function(x, d=NULL, t=NULL, o=FALSE)
{
    if(!is.matrix(x))
        dim(x) <- c(length(x), 1L)
    d <- d %||% 1
    t <- t %||% seq.int(max(d))

    ## get polynomial
    x <- poly(x, degree=max(d), raw=!o, simple=TRUE)

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

pcs <- function(x) x %*% svd(x)$v

rkn <- function(X)  scale(apply(X, 2, rank))
