#'
#' evalcovar.R
#'
#'   evaluate covariate values at data points and at pixels
#'
#' $Revision: 1.23 $ $Date: 2017/08/17 09:20:00 $
#'

evalCovar <- function(model, covariate, ...) {
  UseMethod("evalCovar")
}

evalCovar.ppm <- local({

  evalCovar.ppm <- function(model, covariate, ...,
                            lambdatype=c("cif", "trend", "intensity"),
                            dimyx=NULL, eps=NULL,
                            interpolate=TRUE, jitter=TRUE, 
                            modelname=NULL, covname=NULL,
                            dataname=NULL, subset=NULL) {
    lambdatype <- match.arg(lambdatype)
    #' evaluate covariate values at data points and at pixels
    csr <- is.poisson.ppm(model) && is.stationary.ppm(model)
    #' determine names
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

    #' explicit control of pixel resolution
    if(!is.null(dimyx) || !is.null(eps))
      W <- as.mask(W, dimyx=dimyx, eps=eps)

    if(!is.null(subset)) {
      #' restrict to subset if required
      X <- X[subset]
      W <- W[subset, drop=FALSE]
    }
    
    #' evaluate covariate 
    if(is.character(covariate)) {
      #' One of the characters 'x' or 'y'
      #' Turn it into a function.
      ns <- length(covariate)
      if(ns == 0) stop("covariate is empty")
      if(ns > 1) stop("more than one covariate specified")
      covname <- covariate
      covariate <- switch(covariate,
                          x=xcoordfun,
                          y=ycoordfun,
                          stop(paste("Unrecognised covariate",
                                     dQuote(covariate))))
    } 
  
    if(!is.marked(model)) {
      #' ...................  unmarked .......................
      if(is.im(covariate)) {
        type <- "im"
        if(!interpolate) {
          #' look up covariate values 
          ZX <- safelookup(covariate, X)
        } else {
          #' evaluate at data points by interpolation
          ZX <- interp.im(covariate, X$x, X$y)
          #' fix boundary glitches
          if(any(uhoh <- is.na(ZX)))
            ZX[uhoh] <- safelookup(covariate, X[uhoh])
        }
        #' covariate values for pixels inside window
        Z <- covariate[W, drop=FALSE]
        #' corresponding mask
        W <- as.owin(Z)
      } else if(is.function(covariate)) {
        type <- "function"
        #' evaluate exactly at data points
        ZX <- covariate(X$x, X$y)
        if(!all(is.finite(ZX)))
          warning("covariate function returned NA or Inf values")
        #' window
        W <- as.mask(W)
        #' covariate in window
        Z <- as.im(covariate, W=W)
        #' collapse function body to single string
        covname <- singlestring(covname)
      } else if(is.null(covariate)) {
        stop("The covariate is NULL", call.=FALSE)
      } else stop(paste("The covariate should be",
                        "an image, a function(x,y)",
                        "or one of the characters",
                        sQuote("x"), "or", sQuote("y")),
                  call.=FALSE)
      #' values of covariate in window
      Zvalues <- as.vector(Z[W, drop=TRUE])
      #' corresponding fitted [conditional] intensity values
      lambda <- as.vector(predict(model, locations=W,
                                  type=lambdatype)[W, drop=TRUE])
      #' pixel area (constant)
      pixelarea <- with(Z, xstep * ystep)
    } else {
      #' ...................  marked .......................
      if(!is.multitype(model))
        stop("Only implemented for multitype models (factor marks)")
      marx <- marks(X, dfok=FALSE)
      possmarks <- levels(marx)
      npts <- npoints(X)
      #' single image: replicate 
      if(is.im(covariate)) {
        covariate <- rep(list(covariate), times=length(possmarks))
        names(covariate) <- as.character(possmarks)
      }
      #'
      if(is.list(covariate) && all(unlist(lapply(covariate, is.im)))) {
        #' list of images
        type <- "im"
        if(length(covariate) != length(possmarks))
          stop("Number of images does not match number of possible marks")
        #' evaluate covariate at each data point 
        ZX <- numeric(npts)
        for(k in seq_along(possmarks)) {
          ii <- (marx == possmarks[k])
          covariate.k <- covariate[[k]]
          if(!interpolate) {
            #' look up covariate values 
            values <- safelookup(covariate, X)
          } else {
            #' interpolate
            values <- interp.im(covariate.k, x=X$x[ii], y=X$y[ii])
            #' fix boundary glitches
            if(any(uhoh <- is.na(values)))
              values[uhoh] <- safelookup(covariate.k, X[ii][uhoh])
          }
          ZX[ii] <- values
        }
        #' restrict covariate images to window 
        Z <- lapply(covariate, "[", i=W, drop=FALSE)
        #' extract pixel locations and pixel values
        Zframes <- lapply(Z, as.data.frame)
        #' covariate values at each pixel inside window
        Zvalues <- unlist(lapply(Zframes, getElement, name="value"))
        #' pixel locations 
        locn <- lapply(Zframes, getxy)
        #' tack on mark values
        for(k in seq_along(possmarks))
          locn[[k]] <- cbind(locn[[k]], data.frame(marks=possmarks[k]))
        loc <- do.call(rbind, locn)
        #' corresponding fitted [conditional] intensity values
        lambda <- predict(model, locations=loc, type=lambdatype)
        #' pixel areas
        pixelarea <- rep(sapply(Z, pixarea), sapply(Z, npixdefined))
      } else if(is.function(covariate)) {
        type <- "function"
        #' evaluate exactly at data points
        ZX <- covariate(X$x, X$y, marx)
        #' same window
        W <- as.mask(W)
        #' covariate in window
        Z <- list()
        for(k in seq_along(possmarks))
          Z[[k]] <- as.im(functioncaller, m=possmarks[k], f=covariate, W=W)
        #' functioncaller: function(x,y,m,f) { f(x,y,m) }
        Zvalues <- unlist(lapply(Z, pixelvalues))
        #' corresponding fitted [conditional] intensity values
        lambda <- predict(model, locations=W, type=lambdatype)
        lambda <- unlist(lapply(lambda, pixelvalues))
        if(length(lambda) != length(Zvalues))
          stop("Internal error: length(lambda) != length(Zvalues)")
        #' collapse function body to single string
        covname <- singlestring(covname)
        #' pixel areas
        pixelarea <- rep(sapply(Z, pixarea), sapply(Z, npixdefined))
      } else if(is.null(covariate)) {
        stop("The covariate is NULL", call.=FALSE)
      } else stop(paste("For a multitype point process model,",
                        "the covariate should be an image, a list of images,",
                        "a function(x,y,m)", 
                        "or one of the characters",
                        sQuote("x"), "or", sQuote("y")),
                  call.=FALSE)
    }    
    #' ..........................................................

    #' apply jittering to avoid ties
    if(jitter) {
      nX <- length(ZX)
      dZ <- 0.3 * quantile(diff(sort(unique(c(ZX, Zvalues)))), 1/min(20, nX))
      ZX <- ZX + rnorm(nX, sd=dZ)
      Zvalues <- Zvalues + rnorm(length(Zvalues), sd=dZ)
    }

    lambdaname <- if(is.poisson(model)) "intensity" else lambdatype
    lambdaname <- paste("the fitted", lambdaname)
    check.finite(lambda, xname=lambdaname, usergiven=FALSE)
    check.finite(Zvalues, xname="the covariate", usergiven=TRUE)

    #' lambda values at data points
    lambdaX <- predict(model, locations=X, type=lambdatype)
    
    #' wrap up 
    values <- list(Zimage    = Z,
                   Zvalues   = Zvalues,
                   lambda    = lambda,
                   lambdaX   = lambdaX,
                   weights   = pixelarea,
                   ZX        = ZX,
                   type      = type)
    return(list(values=values, info=info))
  }

  xcoordfun <- function(x,y,m){x}
  ycoordfun <- function(x,y,m){y}

  pixarea <- function(z) { z$xstep * z$ystep }
  npixdefined <- function(z) { sum(!is.na(z$v)) }
  functioncaller <- function(x,y,m,f) { f(x,y,m) }
  pixelvalues <- function(z) { as.data.frame(z)[,3L] }
  getxy <- function(z) { z[,c("x","y")] }
  
  evalCovar.ppm
})

evalCovar.lppm <- local({

  evalCovar.lppm <- function(model, covariate, ...,
                             lambdatype=c("cif", "trend", "intensity"),
                             eps=NULL, nd=1000,
                             interpolate=TRUE, jitter=TRUE, 
                             modelname=NULL, covname=NULL,
                             dataname=NULL, subset=NULL) {
    lambdatype <- match.arg(lambdatype)
    #' evaluate covariate values at data points and at pixels
    csr <- is.poisson(model) && is.stationary(model)

    #' determine names
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

    #' convert character covariate to function
    if(is.character(covariate)) {
      #' One of the characters 'x' or 'y'
      #' Turn it into a function.
      ns <- length(covariate)
      if(ns == 0) stop("covariate is empty")
      if(ns > 1) stop("more than one covariate specified")
      covname <- covariate
      covariate <- switch(covariate,
                          x=xcoordfun,
                          y=ycoordfun,
                          stop(paste("Unrecognised covariate",
                                     dQuote(covariate))))
    }
  
    #' extract model components
    X <- model$X
    fit <- model$fit
    #'
    L <- as.linnet(X)
    Q <- quad.ppm(fit)
    #' restrict to subset if required
    if(!is.null(subset)) {
      X <- X[subset]
      Q <- Q[subset]
    }
    isdat <- is.data(Q)
    U <- union.quad(Q)
    wt <- w.quad(Q)
  
    #' evaluate covariate
    if(!is.marked(model)) {
      #' ...................  unmarked .......................
      if(is.im(covariate)) {
        if(inherits(covariate, "linim")) {
          type <- "linim"
          Zimage <- covariate
        } else {
          type <- "im"
          Zimage <- as.linim(covariate, L)
        }
        if(!interpolate) {
          #' look up covariate values at quadrature points
          Zvalues <- safelookup(covariate, U)
        } else {
          #' evaluate at quadrature points by interpolation
          Zvalues <- interp.im(covariate, U$x, U$y)
          #' fix boundary glitches
          if(any(uhoh <- is.na(Zvalues)))
            Zvalues[uhoh] <- safelookup(covariate, U[uhoh])
        }
        #' extract data values
        ZX <- Zvalues[isdat]
      } else if(is.function(covariate)) {
        type <- "function"
        Zimage <- as.linim(covariate, L)
        #' evaluate exactly at quadrature points
        Zvalues <- covariate(U$x, U$y)
        if(!all(is.finite(Zvalues)))
          warning("covariate function returned NA or Inf values")
        #' extract data values
        ZX <- Zvalues[isdat]
        #' collapse function body to single string
        covname <- singlestring(covname)
      } else if(is.null(covariate)) {
        stop("The covariate is NULL", call.=FALSE)
      } else stop(paste("The covariate should be",
                        "an image, a function(x,y)",
                        "or one of the characters",
                        sQuote("x"), "or", sQuote("y")),
                  call.=FALSE)
      #' corresponding fitted [conditional] intensity values
      lambda <- as.vector(predict(model, locations=U, type=lambdatype))
    } else {
      #' ...................  marked .......................
      if(!is.multitype(model))
      stop("Only implemented for multitype models (factor marks)")
      marx <- marks(U, dfok=FALSE)
      possmarks <- levels(marx)
      #' single image: replicate 
      if(is.im(covariate)) {
        covariate <- rep(list(covariate), length(possmarks))
        names(covariate) <- possmarks
      }
      #'
      if(is.list(covariate) && all(unlist(lapply(covariate, is.im)))) {
        #' list of images
        if(length(covariate) != length(possmarks))
          stop("Number of images does not match number of possible marks")
        #' determine type of data
        islinim <- unlist(lapply(covariate, inherits, what="linim"))
        type <- if(all(islinim)) "linim" else "im"
        Zimage <- covariate
        Zimage[!islinim] <- lapply(Zimage[!islinim], as.linim, L=L)
        #' evaluate covariate at each data point by interpolation
        Zvalues <- numeric(npoints(U))
        for(k in seq_along(possmarks)) {
          ii <- (marx == possmarks[k])
          covariate.k <- covariate[[k]]
          if(!interpolate) {
            #' direct lookup
            values <- safelookup(covariate.k, U[ii])
          } else {
            #' interpolation
            values <- interp.im(covariate.k, x=U$x[ii], y=U$y[ii])
            #' fix boundary glitches
            if(any(uhoh <- is.na(values)))
              values[uhoh] <- safelookup(covariate.k, U[ii][uhoh])
          }
          Zvalues[ii] <- values
        }
        #' extract data values
        ZX <- Zvalues[isdat]
        #' corresponding fitted [conditional] intensity values
        lambda <- predict(model, locations=U, type=lambdatype)
        if(length(lambda) != length(Zvalues))
          stop("Internal error: length(lambda) != length(Zvalues)")
      } else if(is.function(covariate)) {
        type <- "function"
        #' evaluate exactly at quadrature points
        Zvalues <- covariate(U$x, U$y, marx)
        #' extract data values
        ZX <- Zvalues[isdat]
        #' corresponding fitted [conditional] intensity values
        lambda <- predict(model, locations=U, type=lambdatype)
        if(length(lambda) != length(Zvalues))
          stop("Internal error: length(lambda) != length(Zvalues)")
        #' images
        Zimage <- list()
        for(k in seq_along(possmarks))
          Zimage[[k]] <- as.linim(functioncaller, L=L, m=possmarks[k],
                                  f=covariate)
        #' collapse function body to single string
        covname <- singlestring(covname)
      } else if(is.null(covariate)) {
        stop("The covariate is NULL", call.=FALSE)
      } else stop(paste("For a multitype point process model,",
                        "the covariate should be an image, a list of images,",
                        "a function(x,y,m)", 
                        "or one of the characters",
                        sQuote("x"), "or", sQuote("y")),
                  call.=FALSE)
    }    
    #' ..........................................................

    #' apply jittering to avoid ties
    if(jitter) {
      nX <- length(ZX)
      dZ <- 0.3 * quantile(diff(sort(unique(c(ZX, Zvalues)))), 1/min(20, nX))
      ZX <- ZX + rnorm(nX, sd=dZ)
      Zvalues <- Zvalues + rnorm(length(Zvalues), sd=dZ)
    }

    lambdaname <- if(is.poisson(model)) "intensity" else lambdatype
    lambdaname <- paste("the fitted", lambdaname)
    check.finite(lambda, xname=lambdaname, usergiven=FALSE)
    check.finite(Zvalues, xname="the covariate", usergiven=TRUE)

    #' restrict image to subset 
    if(!is.null(subset))
      Zimage <- Zimage[subset, drop=FALSE]

    #' lambda values at data points
    lambdaX <- predict(model, locations=X, type=lambdatype)
    
    #' wrap up 
    values <- list(Zimage    = Zimage,
                   Zvalues   = Zvalues,
                   lambda    = lambda,
                   lambdaX   = lambdaX,
                   weights   = wt,
                   ZX        = ZX,
                   type      = type)
    return(list(values=values, info=info))
  }

  xcoordfun <- function(x,y,m){x}
  ycoordfun <- function(x,y,m){y}
  functioncaller <- function(x,y,m,f) { f(x,y,m) }

  evalCovar.lppm
})
