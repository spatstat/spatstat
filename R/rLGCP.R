#
#   rLGCP.R
#
#   simulation of log-Gaussian Cox process
#
#   original code by Abdollah Jalilian
#
#  $Revision: 1.17 $    $Date: 2015/07/18 03:12:53 $
#

rLGCP <- local({

  rLGCP <- function(model="exp", mu = 0, param = NULL, ...,
                    win=NULL, saveLambda=TRUE, nsim=1, drop=TRUE) {
    ## validate
    if (!(is.numeric(mu) || is.function(mu) || is.im(mu))) 
      stop(paste(sQuote("mu"), "must be a constant, a function or an image"))
    if (is.numeric(mu) && !(length(mu) == 1)) 
      stop(paste(sQuote("mu"), "must be a single number"))
    ## check for outdated usage
    if(!all(nzchar(names(param))))
      stop("Outdated syntax of argument 'param' to rLGCP", call.=FALSE)
    ## 
    do.rLGCP(model=model, mu=mu, param=param, ...,
             win=win, saveLambda=saveLambda, nsim=nsim, drop=drop)
  }

  do.rLGCP <- function(model="exp", mu = 0, param = NULL, ...,
                       win=NULL, saveLambda=TRUE,
                       eps = NULL, dimyx = NULL, xy = NULL,
                       modelonly=FALSE, nsim=1, drop=TRUE) {
    ## make RF model object from RandomFields package
    kraever("RandomFields")
    ## get the 'model generator'
    modelname <- if(model == "exponential") "exp" else model
    modgen <- try(getExportedValue("RandomFields", 
                                   paste0("RM", modelname)),
                  silent=TRUE)
    if(inherits(modgen, "try-error") ||
       !inherits(modgen, "RMmodelgenerator"))
      stop(paste("Model", sQuote(model), "is not recognised"))
    ## now create a RandomFields 'model' object
    rfmodel <- do.call(modgen, append(as.list(param), list(...)))

    ## secret exit
    if(modelonly)
      return(rfmodel)

    ## simulation window
    win.given <- !is.null(win)
    mu.image <- is.im(mu)
    win <- if(win.given) as.owin(win) else if(mu.image) as.owin(mu) else owin()
  
    if(win.given && mu.image && !is.subset.owin(win, as.owin(mu)))
      stop(paste("The spatial domain of the pixel image", sQuote("mu"),
                 "does not cover the simulation window", sQuote("win")))

    ## convert win to a mask
    w <- as.mask(w=win, eps=eps, dimyx=dimyx, xy=xy)
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

    stopifnot(nsim >= 1)
    result <- vector(mode="list", length=nsim)
    for(i in 1:nsim) {
      ## generate zero-mean Gaussian random field
      spc <- RandomFields::RFoptions()$general$spConform
      if(spc) RandomFields::RFoptions(spConform=FALSE)
      z <- RandomFields::RFsimulate(rfmodel, xcol, yrow, grid = TRUE)
      if(spc) RandomFields::RFoptions(spConform=TRUE)

      ## convert to log-Gaussian image
      logLambda <- muxy + z
      Lambda <- matrix(exp(logLambda), nrow=dim[1], ncol=dim[2], byrow=TRUE)
      Lambda <- as.im(Lambda, W=w)
      ## generate Poisson points
      X <- rpoispp(Lambda)[win]
      ## 
      if(saveLambda)
        attr(X, "Lambda") <- Lambda
      result[[i]] <- X
    }
    if(drop && nsim == 1)
      return(result[[1]])
    names(result) <- paste("Simulation", 1:nsim)
    return(as.solist(result))
  }

  rLGCP
})
