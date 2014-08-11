#
# evalcovar.R
#
#   evaluate covariate values at data points and at pixels
#
# $Revision: 1.9 $ $Date: 2013/01/31 07:55:17 $
#

evalCovar <- function(model, covariate, ...) {
  UseMethod("evalCovar")
}

evalCovar.ppm <- function(model, covariate, ...,
                          dimyx=NULL, eps=NULL,
                          jitter=TRUE, 
                          modelname=NULL, covname=NULL,
                          dataname=NULL) {
  # evaluate covariate values at data points and at pixels
  csr <- is.poisson.ppm(model) && is.stationary.ppm(model)

  # determine names
  if(is.null(modelname))
    modelname <- if(csr) "CSR" else short.deparse(substitute(model))
  if(is.null(covname)) {
    covname <- singlestring(short.deparse(substitute(covariate)))
    if(is.character(covariate)) covname <- covariate
  }
  if(is.null(dataname))
    dataname <- model$Qname

  info <-  list(modelname=modelname, covname=covname,
                dataname=dataname, csr=csr,
                spacename="two dimensions")
  
  X <- data.ppm(model)
  W <- as.owin(model)

  # explicit control of pixel resolution
  if(!is.null(dimyx) || !is.null(eps))
    W <- as.mask(W, dimyx=dimyx, eps=eps)

  # evaluate covariate 
  if(is.character(covariate)) {
    # One of the characters 'x' or 'y'
    # Turn it into a function.
    ns <- length(covariate)
    if(ns == 0) stop("covariate is empty")
    if(ns > 1) stop("more than one covariate specified")
    covname <- covariate
    covariate <- switch(covariate,
                     x=function(x,y,m){x},
                     y=function(x,y,m){y},
                     stop(paste("Unrecognised covariate", dQuote(covariate))))
  } 
  
  if(!is.marked(model)) {
    # ...................  unmarked .......................
    if(is.im(covariate)) {
      type <- "im"
      # evaluate at data points by interpolation
      ZX <- interp.im(covariate, X$x, X$y)
      # fix boundary glitches
      if(any(uhoh <- is.na(ZX)))
        ZX[uhoh] <- safelookup(covariate, X[uhoh])
      # covariate values for pixels inside window
      Z <- covariate[W, drop=FALSE]
      # corresponding mask
      W <- as.owin(Z)
    } else if(is.function(covariate)) {
      type <- "function"
      # evaluate exactly at data points
      ZX <- covariate(X$x, X$y)
      if(!all(is.finite(ZX)))
        warning("covariate function returned NA or Inf values")
      # window
      W <- as.mask(W)
      # covariate in window
      Z <- as.im(covariate, W=W)
      # collapse function body to single string
      covname <- singlestring(covname)
    } else stop(paste("The covariate should be",
                      "an image, a function(x,y)",
                      "or one of the characters",
                      sQuote("x"), "or", sQuote("y")))
    # values of covariate in window
    Zvalues <- as.vector(Z[W, drop=TRUE])
    # corresponding fitted intensity values
    lambda <- as.vector(predict(model, locations=W)[W, drop=TRUE])
    # pixel area (constant)
    pixelarea <- with(Z, xstep * ystep)
  } else {
    # ...................  marked .......................
    if(!is.multitype(model))
      stop("Only implemented for multitype models (factor marks)")
    marx <- marks(X, dfok=FALSE)
    possmarks <- levels(marx)
    npts <- npoints(X)
    # single image: replicate 
    if(is.im(covariate))
      covariate <- lapply(possmarks, function(x,v){v}, v=covariate)
    #
    if(is.list(covariate) && all(unlist(lapply(covariate, is.im)))) {
      # list of images
      type <- "im"
      if(length(covariate) != length(possmarks))
        stop("Number of images does not match number of possible marks")
      # evaluate covariate at each data point by interpolation
      ZX <- numeric(npts)
      for(k in seq_along(possmarks)) {
        ii <- (marx == possmarks[k])
        covariate.k <- covariate[[k]]
        values <- interp.im(covariate.k, x=X$x[ii], y=X$y[ii])
        # fix boundary glitches
        if(any(uhoh <- is.na(values)))
          values[uhoh] <- safelookup(covariate.k, X[ii][uhoh])
        ZX[ii] <- values
      }
      # restrict covariate images to window 
      Z <- lapply(covariate, function(x,W){x[W, drop=FALSE]}, W=W)
      # extract pixel locations and pixel values
      Zframes <- lapply(Z, as.data.frame)
      # covariate values at each pixel inside window
      Zvalues <- unlist(lapply(Zframes, function(df) { df[ , 3] }))
      # pixel locations 
      locn <- lapply(Zframes, function(df) { df[ , 1:2] })
      # tack on mark values
      for(k in seq_along(possmarks))
        locn[[k]] <- cbind(locn[[k]], data.frame(marks=possmarks[k]))
      loc <- do.call("rbind", locn)
      # corresponding fitted intensity values
      lambda <- predict(model, locations=loc)
      # pixel areas
      pixelarea <- unlist(lapply(Z, function(z) {
        with(z, rep(xstep * ystep, sum(!is.na(v))))
      }))
    } else if(is.function(covariate)) {
      type <- "function"
      # evaluate exactly at data points
      ZX <- covariate(X$x, X$y, marx)
      # same window
      W <- as.mask(W)
      # covariate in window
      Z <- list()
      g <- function(x,y,m,f) { f(x,y,m) }
      for(k in seq_along(possmarks))
        Z[[k]] <- as.im(g, m=possmarks[k], f=covariate, W=W)
      Zvalues <- unlist(lapply(Z, function(z) { as.data.frame(z)[,3] }))
      # corresponding fitted intensity values
      lambda <- predict(model, locations=W)
      lambda <- unlist(lapply(lambda, function(z) { as.data.frame(z)[,3] }))
      if(length(lambda) != length(Zvalues))
        stop("Internal error: length(lambda) != length(Zvalues)")
      # collapse function body to single string
      covname <- singlestring(covname)
      # pixel areas
      pixelarea <- unlist(lapply(Z, function(z) {
        with(z, rep(xstep * ystep, sum(!is.na(v))))
      }))
    } else stop(paste("For a multitype point process model,",
                      "the covariate should be an image, a list of images,",
                      "a function(x,y,m)", 
                      "or one of the characters",
                      sQuote("x"), "or", sQuote("y")))
  }    
  # ..........................................................

  # apply jittering to avoid ties
  if(jitter) {
    nX <- length(ZX)
    dZ <- 0.3 * quantile(diff(sort(unique(c(ZX, Zvalues)))), 1/min(20, nX))
    ZX <- ZX + rnorm(nX, sd=dZ)
    Zvalues <- Zvalues + rnorm(length(Zvalues), sd=dZ)
  }

  check.finite(lambda, xname="the fitted intensity lambda", usergiven=FALSE)
  check.finite(Zvalues, xname="the covariate", usergiven=TRUE)
  
  # wrap up 
  values <- list(Zimage    = Z,
                 Zvalues   = Zvalues,
                 lambda    = lambda,
                 weights   = pixelarea,
                 ZX        = ZX,
                 type      = type)
  return(list(values=values, info=info))
}

evalCovar.lppm <- function(model, covariate, ...,
                           eps=NULL, nd=1000,
                           jitter=TRUE, 
                           modelname=NULL, covname=NULL,
                           dataname=NULL) {
  # evaluate covariate values at data points and at pixels
  csr <- is.poisson(model) && is.stationary(model)

  # determine names
  if(is.null(modelname))
    modelname <- if(csr) "CSR" else short.deparse(substitute(model))
  if(is.null(covname)) {
    covname <- singlestring(short.deparse(substitute(covariate)))
    if(is.character(covariate)) covname <- covariate
  }
  if(is.null(dataname))
    dataname <- model$Xname
  info <-  list(modelname=modelname, covname=covname,
                dataname=dataname, csr=csr,
                spacename="linear network")

  # convert character covariate to function
  if(is.character(covariate)) {
    # One of the characters 'x' or 'y'
    # Turn it into a function.
    ns <- length(covariate)
    if(ns == 0) stop("covariate is empty")
    if(ns > 1) stop("more than one covariate specified")
    covname <- covariate
    covariate <- switch(covariate,
                     x=function(x,y,m){x},
                     y=function(x,y,m){y},
                     stop(paste("Unrecognised covariate", dQuote(covariate))))
  }
  
  # extract model components
  X <- model$X
  fit <- model$fit
  #
  L <- as.linnet(X)
  Q <- quad.ppm(fit)
  isdat <- is.data(Q)
  U <- union.quad(Q)
  wt <- w.quad(Q)
  
  # evaluate covariate
  if(!is.marked(model)) {
    # ...................  unmarked .......................
    if(is.im(covariate)) {
      if(inherits(covariate, "linim")) {
        type <- "linim"
        Zimage <- covariate
      } else {
        type <- "im"
        Zimage <- as.linim(covariate, L)
      }
      # evaluate at quadrature points by interpolation
      Zvalues <- interp.im(covariate, U$x, U$y)
      # fix boundary glitches
      if(any(uhoh <- is.na(Zvalues)))
        Zvalues[uhoh] <- safelookup(covariate, U[uhoh])
      # extract data values
      ZX <- Zvalues[isdat]
    } else if(is.function(covariate)) {
      type <- "function"
      Zimage <- as.linim(covariate, L)
      # evaluate exactly at quadrature points
      Zvalues <- covariate(U$x, U$y)
      if(!all(is.finite(Zvalues)))
        warning("covariate function returned NA or Inf values")
      # extract data values
      ZX <- Zvalues[isdat]
      # collapse function body to single string
      covname <- singlestring(covname)
    } else stop(paste("The covariate should be",
                      "an image, a function(x,y)",
                      "or one of the characters",
                      sQuote("x"), "or", sQuote("y")))
    # corresponding fitted intensity values
    lambda <- as.vector(predict(model, locations=U))
  } else {
    # ...................  marked .......................
    if(!is.multitype(model))
      stop("Only implemented for multitype models (factor marks)")
    marx <- marks(U, dfok=FALSE)
    possmarks <- levels(marx)
    npts <- npoints(X)
    # single image: replicate 
    if(is.im(covariate))
      covariate <- lapply(possmarks, function(x,v){v}, v=covariate)
    #
    if(is.list(covariate) && all(unlist(lapply(covariate, is.im)))) {
      # list of images
      if(length(covariate) != length(possmarks))
        stop("Number of images does not match number of possible marks")
      # determine type of data
      islinim <- unlist(lapply(covariate, inherits, what="linim"))
      type <- if(all(islinim)) "linim" else "im"
      Zimage <- covariate
      Zimage[!islinim] <- lapply(Zimage[!islinim], as.linim, L=L)
      # evaluate covariate at each data point by interpolation
      Zvalues <- numeric(npoints(U))
      for(k in seq_along(possmarks)) {
        ii <- (marx == possmarks[k])
        covariate.k <- covariate[[k]]
        values <- interp.im(covariate.k, x=U$x[ii], y=U$y[ii])
        # fix boundary glitches
        if(any(uhoh <- is.na(values)))
          values[uhoh] <- safelookup(covariate.k, U[ii][uhoh])
        Zvalues[ii] <- values
      }
      # extract data values
      ZX <- Zvalues[isdat]
      # corresponding fitted intensity values
      lambda <- predict(model, locations=U)
      if(length(lambda) != length(Zvalues))
        stop("Internal error: length(lambda) != length(Zvalues)")
    } else if(is.function(covariate)) {
      type <- "function"
      # evaluate exactly at quadrature points
      Zvalues <- covariate(U$x, U$y, marx)
      # extract data values
      ZX <- Zvalues[isdat]
      # corresponding fitted intensity values
      lambda <- predict(model, locations=U)
      if(length(lambda) != length(Zvalues))
        stop("Internal error: length(lambda) != length(Zvalues)")
      # images
      Zimage <- list()
      g <- function(x,y,m,f) { f(x,y,m) }
      for(k in seq_along(possmarks))
        Zimage[[k]] <- as.linim(g, L=L, m=possmarks[k], f=covariate)
      # collapse function body to single string
      covname <- singlestring(covname)
    } else stop(paste("For a multitype point process model,",
                      "the covariate should be an image, a list of images,",
                      "a function(x,y,m)", 
                      "or one of the characters",
                      sQuote("x"), "or", sQuote("y")))
  }    
  # ..........................................................

  # apply jittering to avoid ties
  if(jitter) {
    nX <- length(ZX)
    dZ <- 0.3 * quantile(diff(sort(unique(c(ZX, Zvalues)))), 1/min(20, nX))
    ZX <- ZX + rnorm(nX, sd=dZ)
    Zvalues <- Zvalues + rnorm(length(Zvalues), sd=dZ)
  }

  check.finite(lambda, xname="the fitted intensity lambda", usergiven=FALSE)
  check.finite(Zvalues, xname="the covariate", usergiven=TRUE)

  # wrap up 
  values <- list(Zimage    = Zimage,
                 Zvalues   = Zvalues,
                 lambda    = lambda,
                 weights   = wt,
                 ZX        = ZX,
                 type      = type)
  return(list(values=values, info=info))
}

