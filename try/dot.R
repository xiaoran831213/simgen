#' Schur complement
#'
#' Given a full matrix \code{X} and a target block \code{C} within,
#' 
#' X = [A  B]
#'     [B' C],
#'
#' calculate the Schur complement \eqn{X/C = A - BC^{-1}B'}
#' 
#' @param X the whole matrix
#' @param C mask or index of the target block
#' @return the Schur complement of block C of matrix X
#'
#' @examples
#' ## get genotype and covariate matrices
#' gno <- readRDS(system.file("extdata", 'rs208294_gno.rds', package="dotgen"))
#' cvr <- readRDS(system.file("extdata", 'rs208294_cvr.rds', package="dotgen"))
#' 
#' ## overall correlation is a 56 x 56 positive definite matrix
#' X <- cor(cbind(gno, cvr))
#'
#' ## Schur complement for the covariate block
#' C <- 1:ncol(cvr) + ncol(gno)      # block of covariate
#' A <- 1:ncol(gno)                  # block of genotype
#' res1 <- solve(solve(X)[A, A])     # method 1
#' res2 <- scp(X, C)                 # method 2
#'
#' ## check equality
#' stopifnot(all.equal(res1, res2))  # must be TRUE
#' 
#' @noRd
scp <- function(X, C)
{
    ## C can be names, numbers, or Boolean masks, turn them into index
    if(is.character(C))
        C <- colnames(X) %in% C
    if(is.logical(C))
        C <- which(C)
    i <- unname(C)
    
    ## A <- X[-C, -C]
    ## B <- X[-C, +C]
    ## B'<- X[+C, -C]
    ## C <- X[+C, +C]

    ## A - B C^{-1} B'
    X[-i, -i] - X[-i, +i] %*% solve(X[+i, +i], X[+i, -i])
}

#' Multivariate negative square root of a PD matrix
#'
#' Decorrelate assocation test statistics between multiple phenotypes and g-variants.
#' 
#' @param C M x M matrix of correlation among M variants adjusted for covariants.
#' @param D Q x Q matrix of correlation among Q multivariate.
nsp <- function(C, D=NULL, tol.egv=NULL, ...)
{
    tol.egv <- tol.egv %||%  sqrt(.Machine$double.eps)
    D <- D %||% 1

    C <- eigen(C, TRUE)
    D <- eigen(D, TRUE)
    
    d1 <- D$values
    d2 <- C$values
    i1 <- d1 > max(d1) * tol.egv
    i2 <- d2 > max(d2) * tol.egv
    d1 <- d1[i1]
    d2 <- d2[i2]
    u1 <- D$vectors[, i1]
    u2 <- C$vectors[, i2]
    d <- kronecker(d1, d2)
    u <- kronecker(u1, u2)
    
    ## d <- kronecker(D$values, C$values)
    ## u <- kronecker(D$vectors, C$vectors)
    
    ## positive eigen values
    ## . <- d > max(d) * tol.egv
    ## if(!all(.))
    ## {
    ##     d <- d[  .]
    ##     u <- u[, .]
    ## }
    L <- length(d)              # effective number of eigen
    
    ## square root
    dim(d) <- NULL
    d <- sqrt(1/d)
    H <- u %*% (d * t(u))       # U diag(d) U'
    H <- 0.5 * (H + t(H))
    
    list(H=H, L=L)
}

#' Multivariate Decorrelation by Orthogonal Transformation
#'
#' Decorrelate assocation test statistics between multiple phenotypes and g-variants.
#' 
#' @param Z Q x M matrix of z-scores between Q traits and M variants.
#' @param C M x M matrix of correlation among M variants adjusted for covariants.
#' @param D Q x Q matrix of correlation among Q traits.
#' @param tol.cor tolerance threshold for the largest correlation absolute value.
#' @param tol.egv tolerance threshold for the smallest eigenvalue.
#' @param ... additional parameters.
mdt <- function(Z, C, D=NULL, tol.cor=NULL, tol.egv=NULL, use.egv=0, ...)
{
    if(is.null(tol.cor))
        tol.cor <- sqrt(.Machine$double.eps)
    if(is.null(tol.egv))
        tol.egv <- sqrt(.Machine$double.eps)
    
    ## trim collinear variants
    m <- dvt(C, tol.cor)
    M <- sum(m)                      # effective number of variants
    Z <- c(Z[m, ])
    C <- C[m, m]
    
    ret <- nsp(C, D, tol.egv=tol.egv, ...)
    H <- ret$H
    L <- ret$L
    
    X <- H %*% c(Z)

    P <- 1 - pchisq(sum(X^2), L)
    c(list(X=X, P=P), ret, M=M)
}

dvt <- function(C, tol.cor=NULL, ...)
{
    if(is.null(tol.cor))
        tol.cor <- sqrt(.Machine$double.eps)
    
    C <- abs(C)
    C[upper.tri(C, TRUE)] <- 0
    apply(C, 2, function(.) all(. < 1 - tol.cor))
}
