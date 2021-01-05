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
