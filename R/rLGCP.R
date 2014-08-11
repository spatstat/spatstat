#
#   rLGCP.R
#
#   simulation of log-Gaussian Cox process
#
#   original code by Abdollah Jalilian
#
#  $Revision: 1.7 $    $Date: 2013/01/13 04:01:55 $
#

rLGCP <-
  function(model="exponential", mu = 0, param = NULL, ..., win=NULL)
{
  if(!missing(mu)) {
    if (!(is.numeric(mu) || is.function(mu) || is.im(mu))) 
      stop(paste(sQuote("mu"), "must be a constant, a function or an image"))
    if (is.numeric(mu) && !(length(mu) == 1)) 
      stop(paste(sQuote("mu"), "must be a single number"))
  }
  if(!require(RandomFields))
    stop("Simulation of log-Gaussian Cox process requires the package RandomFields")
  win.given <- !is.null(win)
  mu.image <- is.im(mu)
  win <- if(win.given) as.owin(win) else if(mu.image) as.owin(mu) else owin()
  
  if(win.given && mu.image && !is.subset.owin(win, as.owin(mu)))
    stop(paste("The spatial domain of the pixel image", sQuote("mu"),
               "does not cover the simulation window", sQuote("win")))
  
  w <- as.mask(win)
  x <- w$xcol
  y <- w$yrow
  dim <- w$dim
  xy <- expand.grid(x=x, y=y)
  xx <- xy$x
  yy <- xy$y

  muxy <- if(is.numeric(mu)) mu else
          if (is.function(mu)) mu(xx,yy) else
          lookup.im(mu, xx, yy, naok=TRUE, strict=TRUE)
  muxy[is.na(muxy)] <- -Inf

  # generate Gaussian random field
  xgrid <- c(x[1], x[length(x)], w$xstep)
  ygrid <- c(y[1], y[length(y)], w$ystep)
  z <- RandomFields::GaussRF(xgrid, ygrid, grid = TRUE, gridtriple=TRUE,
               model = model, param = param, ...)
  logLambda <- muxy + z
  # convert to log-Gaussian image
  Lambda <- matrix(exp(logLambda), nrow=dim[1], ncol=dim[2], byrow=TRUE)
  Lambda <- as.im(Lambda, W=w)
  # generate Poisson points
  X <- rpoispp(Lambda)[win]
  # 
  attr(X, "Lambda") <- Lambda
  return(X)
}

