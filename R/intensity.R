#
# intensity.R
#
# Code related to intensity and intensity approximations
#
#  $Revision: 1.16 $ $Date: 2016/10/09 02:13:53 $
#

intensity <- function(X, ...) {
  UseMethod("intensity")
}

intensity.ppp <- function(X, ..., weights=NULL) {
  n <- npoints(X)
  a <- area(Window(X))
  if(is.null(weights)) {
    ## unweighted case - for efficiency
    if(is.multitype(X)) {
      mks <- marks(X)
      answer <- as.vector(table(mks))/a
      names(answer) <- levels(mks)
    } else answer <- n/a
    return(answer)
  }
  ## weighted case 
  if(is.numeric(weights)) {
    check.nvector(weights, n)
  } else if(is.expression(weights)) {
    # evaluate expression in data frame of coordinates and marks
    df <- as.data.frame(X)
    pf <- parent.frame()
    eval.weights <- try(eval(weights, envir=df, enclos=pf))
    if(inherits(eval.weights, "try-error"))
      stop("Unable to evaluate expression for weights", call.=FALSE)
    if(!check.nvector(eval.weights, n, fatal=FALSE, warn=TRUE))
      stop("Result of evaluating the expression for weights has wrong format")
    weights <- eval.weights
  } else stop("Unrecognised format for argument 'weights'")
  ##
  if(is.multitype(X)) {
    mks <- marks(X)
    answer <- as.vector(tapply(weights, mks, sum))/a
    answer[is.na(answer)] <- 0
    names(answer) <- levels(mks)
  } else {
    answer <- sum(weights)/a
  }
  return(answer)
}

intensity.splitppp <- function(X, ..., weights=NULL) {
  if(is.null(weights))
    return(sapply(X, intensity.ppp))
  if(is.expression(weights))
    return(sapply(X, intensity.ppp, weights=weights))
  if(is.numeric(weights)) {
    fsplit <- attr(X, "fsplit")
    n <- length(fsplit)
    check.nvector(weights, n)
    result <- mapply(intensity.ppp, X, weights=split(weights, fsplit))
    result <- simplify2array(result, higher=FALSE)
    return(result)
  }
  stop("Unrecognised format for weights")
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
      lam <- sX$trend$value
      if(sX$multitype && sX$no.trend) {
        ## trend is ~1; lam should be replicated for each mark
        lev <- levels(marks(data.ppm(X)))
        lam <- rep(lam, length(lev))
        names(lam) <- lev
      }
      return(lam)
    }
    # Nonstationary Poisson
    return(predict(X, ...))
  }
  # Gibbs process
  if(is.multitype(X))
    stop("Not yet implemented for multitype Gibbs processes")
  # Compute first order term
  if(is.stationary(X)) {
    ## activity parameter
    sX <- summary(X, quick="no variances")
    beta <- sX$trend$value
  } else {
    ## activity function (or values of it, depending on '...')
    beta <- predict(X, ...)
  }
  ## apply approximation
  lambda <- PoisSaddleApp(beta, fitin(X))
  return(lambda)
}

PoisSaddleApp <- function(beta, fi) {
  ## apply Poisson-Saddlepoint approximation
  ## given first order term and fitted interaction
  stopifnot(inherits(fi, "fii"))
  inte <- as.interact(fi)
  if(!identical(inte$family$name, "pairwise"))
    stop("Intensity approximation is only available for pairwise interaction models")
  # Stationary, pairwise interaction
  Mayer <- inte$Mayer
  if(is.null(Mayer))
    stop(paste("Sorry, not yet implemented for", inte$name))
  # interaction coefficients
  co <- with(fi, coefs[Vnames[!IsOffset]])
  # compute second Mayer cluster integral
  G <- Mayer(co, inte)
  if(is.null(G) || !is.finite(G)) 
    stop("Internal error in computing Mayer cluster integral")
  if(G < 0)
    stop(paste("Unable to apply Poisson-saddlepoint approximation:",
               "Mayer cluster integral is negative"))
  ## solve
  if(is.im(beta)) {
    lambda <- if(G == 0) beta else eval.im(LambertW(G * beta)/G)
  } else {
    lambda <- if(G == 0) beta else (LambertW(G * beta)/G)
    if(length(lambda) == 1) lambda <- unname(lambda)
  }
  return(lambda)
}


# Lambert's W-function

LambertW <- local({

  yexpyminusx <- function(y,x){y*exp(y)-x}

  W <- function(x) {
    result <- rep.int(NA_real_, length(x))
    ok <- is.finite(x) & (x >= 0)
    if(requireNamespace("gsl", quietly=TRUE)) {
      result[ok] <- gsl::lambert_W0(x[ok])
    } else {
      for(i in which(ok))
        result[i] <- uniroot(yexpyminusx, c(0, x[i]), x=x[i])$root
    }
    return(result)
  }

  W
})

