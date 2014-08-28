#
#   rLGCP.R
#
#   simulation of log-Gaussian Cox process
#
#   original code by Abdollah Jalilian
#
#  $Revision: 1.10 $    $Date: 2014/08/27 09:50:07 $
#

rLGCP <- function(model="exp", mu = 0, param = NULL, ...,
                  win=NULL, saveLambda=TRUE)
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

  ## convert win to a mask
  dc <- do.call.matched(as.mask, append(list(w=win), list(...)), sieve=TRUE)
  w  <- dc$result
  dotargs <- dc$otherargs
  
  xcol <- w$xcol
  yrow <- w$yrow
  dim <- w$dim
  xy <- expand.grid(x=xcol, y=yrow)
  xx <- xy$x
  yy <- xy$y

  muxy <- if(is.numeric(mu)) mu else
          if (is.function(mu)) mu(xx,yy) else
          lookup.im(mu, xx, yy, naok=TRUE, strict=TRUE)
  muxy[is.na(muxy)] <- -Inf

  ## check for outdated usage
  plist <- as.list(param)
  if(!all(nzchar(names(plist))))
    stop("Outdated syntax of argument 'param' to rLGCP", call.=FALSE)
  
  ## get the 'model generator'
  modelname <- if(model == "exponential") "exp" else model
  modgen <- mget(paste0("RM", modelname), inherits=TRUE,
                 ifnotfound=list(NULL))[[1]]
  if(is.null(modgen) || !inherits(modgen, "RMmodelgenerator"))
    stop(paste("Model", sQuote(modelname), "is not recognised"))
  ## now create a RandomFields 'model' object
  rfmodel <- do.call(modgen, append(plist, dotargs))
  
  ## generate zero-mean Gaussian random field
  spc <- RandomFields::RFoptions()$general$spConform
  if(spc) RandomFields::RFoptions(spConform=FALSE)
  z <- RandomFields::RFsimulate(rfmodel, xcol, yrow, grid = TRUE)
  if(spc) RandomFields::RFoptions(spConform=TRUE)

  ## convert to log-Gaussian image
  logLambda <- muxy + z
  Lambda <- matrix(exp(logLambda), nrow=dim[1], ncol=dim[2], byrow=TRUE)
  Lambda <- as.im(Lambda, W=w)
  # generate Poisson points
  X <- rpoispp(Lambda)[win]
  # 
  if(saveLambda)
    attr(X, "Lambda") <- Lambda
  
  return(X)
}

