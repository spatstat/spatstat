#
#   rLGCP.R
#
#   simulation of log-Gaussian Cox process
#
#   original code by Abdollah Jalilian
#
#  $Revision: 1.21 $    $Date: 2019/07/09 09:52:18 $
#

rLGCP <- local({

  rLGCP <- function(model="exp", mu = 0, param = NULL, ...,
                    win=NULL, saveLambda=TRUE, nsim=1, drop=TRUE) {
    ## validate
    if (is.numeric(mu)) {
      check.1.real(mu, paste("if", sQuote("mu"), "is numeric,"))
    } else if(!is.function(mu) && !is.im(mu))
      stop(paste(sQuote("mu"), "must be a constant, a function or an image"))
    check.1.integer(nsim)
    stopifnot(nsim >= 1)
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
    ## get the 'model generator'
    modgen <- getRandomFieldsModelGen(model)
    ## now create a RandomFields 'model' object
    rfmodel <- do.call(modgen, append(as.list(param), list(...)))
    if(!inherits(rfmodel, "RMmodel"))
      stop("Unable to create RandomFields model object", call.=FALSE)

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
    dimw <- w$dim

    ## evaluate 'mu' at pixels of mask
    if(is.numeric(mu)) {
      muxy <- mu
    } else {
      xy <- rasterxy.mask(w, drop=FALSE)
      xx <- xy$x
      yy <- xy$y
      muxy <- if (is.function(mu)) mu(xx,yy) else
              lookup.im(mu, xx, yy, naok=TRUE, strict=TRUE)
      muxy[is.na(muxy)] <- -Inf
    }
    ## corresponding image template
    Lambda <- as.im(w)

    ## generate 'nsim' realisations of a zero-mean Gaussian random field Z
    spc <- RandomFields::RFoptions()$general$spConform
    if(spc) RandomFields::RFoptions(spConform=FALSE)
    z <- RandomFields::RFsimulate(rfmodel, xcol, yrow, grid = TRUE, n=nsim)
    if(spc) RandomFields::RFoptions(spConform=TRUE)
    ## ensure 3D array
    if(length(dim(z)) == 2) z <- array(z, dim=c(dim(z), 1))

    ## generate realisations of LGCP
    result <- vector(mode="list", length=nsim)
    for(i in 1:nsim) {
      ## Extract i-th realisation of Z; convert to log-Gaussian image
      Lambda$v[] <- exp(muxy + z[,,i])
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
