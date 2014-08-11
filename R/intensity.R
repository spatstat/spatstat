#
# intensity.R
#
# Code related to intensity and intensity approximations
#
#  $Revision: 1.4 $ $Date: 2013/01/13 04:01:07 $
#

intensity <- function(X, ...) {
  UseMethod("intensity")
}

intensity.ppp <- function(X, ...) {
  sX <- summary(X, quick="no variances")
  if(is.multitype(X)) {
    answer <- sX$marks$intensity
    names(answer) <- row.names(sX$marks)
  } else {
    answer <- unname(sX$intensity)
  }
  return(answer)
}

intensity.ppm <- function(X, ...) {
  if(!identical(valid.ppm(X), TRUE)) {
    warning("Model is invalid - projecting it")
    X <- project.ppm(X)
  }
  if(is.poisson(X)) {
    if(is.stationary(X)) {
      # stationary univariate/multivariate Poisson
      sX <- summary(X, quick="no variances")
      return(sX$trend$value)
    }
    # Nonstationary Poisson
    return(predict(X, ...))
  }
  # Gibbs process
  if(!is.stationary(X))
    stop("Not yet implemented for non-stationary Gibbs models")
  if(is.multitype(X))
    stop("Not yet implemented for multitype Gibbs processes")
  inte <- as.interact(X)
  if(!identical(inte$family$name, "pairwise"))
    stop("Intensity approximation is only available for pairwise interaction models")
  # Stationary, pairwise interaction
  Mayer <- inte$Mayer
  if(is.null(Mayer))
    stop(paste("Sorry, not yet implemented for", inte$name))
  # interaction coefficients
  co <- with(fitin(X), coefs[Vnames[!IsOffset]])
  # compute second Mayer cluster integral
  G <- Mayer(co, inte)
  if(is.null(G) || !is.finite(G)) 
    stop("Internal error in computing Mayer cluster integral")
  if(G < 0)
    stop(paste("Unable to apply Poisson-saddlepoint approximation:",
               "Mayer cluster integral is negative"))
  # activity parameter
  sX <- summary(X, quick="no variances")
  beta <- sX$trend$value
  # solve
  lambda <- if(G == 0) rep(0, length(beta)) else LambertW(G * beta)/G
  if(length(lambda) == 1) lambda <- unname(lambda)
  return(lambda)
}

# Lambert's W-function

LambertW <- local({

  yexpyminusx <- function(y,x){y*exp(y)-x}

  W <- function(x) {
    if(require(gsl, quietly=TRUE))
      return(gsl::lambert_W0(x))
    result <- rep(NA_real_, length(x))
    for(i in which(is.finite(x) & (x >= 0)))
      result[i] <- uniroot(yexpyminusx, c(0, x[i]), x=x[i])$root
  return(result)
  }

  W
})

