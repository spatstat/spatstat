#'     bc.R
#' 
#'  Bias correction techniques
#'
#'  $Revision: 1.2 $ $Date: 2016/09/15 02:21:15 $

bc <- function(fit, ...) {
  UseMethod("bc")
}

bc.ppm <- function(fit, ..., nfine=256) {
  stopifnot(is.ppm(fit))
  #
  theta0 <- coef(fit)
  nc <- length(theta0)
  #
  X <- data.ppm(fit)
  Z <- is.data(quad.ppm(fit))
  # evaluate sufficient statistic at data points
  sufX <- model.matrix(fit)[Z, ]
  if(ncol(sufX) != nc)
    stop("Internal error: model.matrix does not match coef(model)")
  # predict on fine grid
  finemask <- as.mask(as.owin(fit), dimyx=nfine)
  lamF <- predict(fit, type="cif", locations=finemask)
  sufF <- model.images(fit, W=finemask)
  if(length(sufF) != nc)
    stop("Internal error: model.images does not match coef(model)")
  # edge correction
  if(fit$correction == "border" && ((rbord <- fit$rbord) > 0)) {
    b <- bdist.pixels(finemask)
    bX <- bdist.points(X)
    excludeU <- eval.im(b < rbord)
    retainX  <- (bX >= rbord)
    sufX <- sufX[retainX, , drop=FALSE]
  } else {
    excludeU <- FALSE
  }
  # compute fine approximation to score
  scoreX <- colSums(sufX)
  scoreW <- numeric(nc)
  for(k in seq_len(nc)) {
    S <- sufF[[k]]
    # infinite values of S may occur and correspond to zero cif
    Slam <- eval.im(ifelse(is.infinite(S) | excludeU, 0, S * lamF))
    scoreW[k] <- integral.im(Slam)
  }
  score <- scoreX - scoreW
  # Newton-Raphson
  Iinv <- vcov(fit, hessian=TRUE)
  theta <- theta0 + Iinv %*% score
  theta <- theta[ , 1, drop=TRUE]
  #
#  return(list(theta0=theta0, theta=theta))
  return(theta)
}

# Richardson extrapolation (generic)

rex <- function(x, r=2, k=1, recursive=FALSE) {
  # x should be a matrix
  # whose columns are successive estimates of a parameter vector
  # obtained using "grid step sizes" t, t/r, t/r^2, ...
  # Estimate from step size t is assumed to converge at rate t^k
  if(!is.matrix(x)) x <- matrix(x, nrow=1)
  if(ncol(x) <= 1) return(x)
  rk <- r^k
  y <- (rk * x[, -1, drop=FALSE] - x[, -ncol(x), drop=FALSE])/(rk - 1)
  if(recursive)
    y <- rex(y, r=r, k=k+1, recursive=TRUE)
  return(y)
}

