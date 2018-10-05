#
#  smooth.ppp.R
#
#  Smooth the marks of a point pattern
# 
#  $Revision: 1.68 $  $Date: 2018/10/05 03:49:46 $
#

# smooth.ppp <- function(X, ..., weights=rep(1, npoints(X)), at="pixels") {
#   .Deprecated("Smooth.ppp", package="spatstat",
#    msg="smooth.ppp is deprecated: use the generic Smooth with a capital S")
#   Smooth(X, ..., weights=weights, at=at)
# }

Smooth <- function(X, ...) {
  UseMethod("Smooth")
}

Smooth.solist <- function(X, ...) {
  solapply(X, Smooth, ...)
}

Smooth.ppp <- function(X, sigma=NULL, ...,
                       weights=rep(1, npoints(X)), at="pixels",
                       adjust=1, varcov=NULL, 
                       edge=TRUE, diggle=FALSE, 
                       kernel="gaussian",
                       scalekernel=is.character(kernel),
                       geometric=FALSE) {
  verifyclass(X, "ppp")
  if(!is.marked(X, dfok=TRUE, na.action="fatal"))
    stop("X should be a marked point pattern", call.=FALSE)
  X <- coerce.marks.numeric(X)
  if(!all(is.finite(as.matrix(marks(X)))))
    stop("Some mark values are Inf, NaN or NA", call.=FALSE)
  at <- pickoption("output location type", at,
                   c(pixels="pixels",
                     points="points"))
  ## trivial case
  if(npoints(X) == 0) {
    cn <- colnames(marks(X))
    nc <- length(cn)
    switch(at,
           points = {
             result <- if(nc == 0) numeric(0) else
             matrix(, 0, nc, dimnames=list(NULL, cn))
           },
           pixels = {
             result <- as.im(NA_real_, Window(X))
             if(nc) {
               result <- as.solist(rep(list(result), nc))
               names(result) <- cn
             }
           })
    return(result)
  }

  ## ensure weights are numeric
  weightsgiven <- !missing(weights) && !is.null(weights) 
  if(weightsgiven) {
    # convert to numeric
    if(is.im(weights)) {
      weights <- safelookup(weights, X) # includes warning if NA
    } else if(is.expression(weights)) 
      weights <- eval(weights, envir=as.data.frame(X), enclos=parent.frame())
    if(length(weights) == 0)
      weightsgiven <- FALSE
  }
  if(weightsgiven) {
    check.nvector(weights, npoints(X))
  } else weights <- NULL

  ## geometric mean smoothing
  if(geometric) 
    return(ExpSmoothLog(X, sigma=sigma, ..., at=at,
                        adjust=adjust, varcov=varcov,
                        kernel=kernel, scalekernel=scalekernel,
                        weights=weights, edge=edge, diggle=diggle))

  ## determine smoothing parameters
  if(scalekernel) {
    ker <- resolve.2D.kernel(sigma=sigma, ...,
                             adjust=adjust, varcov=varcov,
                             kernel=kernel, 
                             x=X, bwfun=bw.smoothppp, allow.zero=TRUE)
    sigma <- ker$sigma
    varcov <- ker$varcov
    adjust <- 1
  } 
  
  ## infinite bandwidth
  if(bandwidth.is.infinite(sigma)) {
    #' uniform estimate
    nX <- npoints(X)
    if(is.null(weights)) weights <- rep(1, nX)
    marx <- marks(X)
    single <- is.null(dim(marx))
    wtmark <- weights * marx 
    totwt <- sum(weights)
    totwtmark <- if(single) sum(wtmark) else colSums(wtmark)
    W <- Window(X)
    switch(at,
           pixels = {
             result <- solapply(totwtmark/totwt, as.im, W=W, ...)
             names(result) <- colnames(marx)
             if(single) result <- result[[1L]]
           },
           points = {
             denominator <- rep(totwt, nX)
             numerator <- rep(totwtmark, each=nX)
             if(!single) numerator <- matrix(numerator, nrow=nX)
             leaveoneout <- resolve.1.default(list(leaveoneout=TRUE), list(...))
             if(leaveoneout) {
               numerator <- numerator - wtmark
               denominator <- denominator - weights
             }
             result <- numerator/denominator
             if(!single)
               colnames(result) <- colnames(marx)
           })
    return(result)
  }

  ## Diggle's edge correction?
  if(diggle && !edge) warning("Option diggle=TRUE overridden by edge=FALSE")
  diggle <- diggle && edge
  ##
  ## cutoff distance (beyond which the kernel value is treated as zero)
  cutoff <- cutoff2Dkernel(kernel, sigma=sigma, varcov=varcov,
                           scalekernel=scalekernel, adjust=adjust, ...,
                           fatal=TRUE)
  ## 
  if(cutoff < minnndist(X)) {
    # very small bandwidth
    leaveoneout <- resolve.1.default("leaveoneout",
                                     list(...), list(leaveoneout=TRUE))
    if(!leaveoneout && at=="points") {
      warning(paste("Bandwidth is close to zero:",
                    "original values returned"))
      Y <- marks(X)
    } else {
      warning(paste("Bandwidth is close to zero:",
                    "nearest-neighbour interpolation performed"))
      Y <- nnmark(X, ..., k=1, at=at)
    }
    return(Y)
  }

  if(diggle) {
    ## absorb Diggle edge correction into weights vector
    edg <- second.moment.calc(X, sigma, what="edge", ...,                                                     varcov=varcov, adjust=adjust,
                              kernel=kernel, scalekernel=scalekernel)
    ei <- safelookup(edg, X, warn=FALSE)
    weights <- if(weightsgiven) weights/ei else 1/ei
    weights[!is.finite(weights)] <- 0
    weightsgiven <- TRUE
  }
  ## rescale weights to avoid numerical gremlins
  if(weightsgiven && ((mw <- median(abs(weights))) > 0))
    weights <- weights/mw

  ## calculate...
  marx <- marks(X)
  if(!is.data.frame(marx)) {
    # ........ vector of marks ...................
    values <- marx
    if(is.factor(values)) {
      warning("Factor valued marks were converted to integers")
      values <- as.numeric(values)
    }
    ## detect constant values
    ra <- range(values, na.rm=TRUE)
    if(diff(ra) == 0) {
      switch(at,
             points = {
               result <- values
             },
             pixels = {
               M <- do.call.matched(as.mask, list(w=as.owin(X), ...))
               result <- as.im(ra[1], M)
             })
    } else {
      switch(at,
             points={
               result <-
                 do.call(smoothpointsEngine,
                         resolve.defaults(list(x=X,
                                               values=values, weights=weights,
                                               sigma=sigma, varcov=varcov,
                                               kernel=kernel,
                                               scalekernel=scalekernel,
                                               edge=FALSE),
                                          list(...)))
             },
             pixels={
               values.weights <- if(weightsgiven) values * weights else values
               numerator <-
                 do.call(density.ppp,
                         resolve.defaults(list(x=X,
                                               at="pixels",
                                               weights = values.weights,
                                               sigma=sigma, varcov=varcov,
                                               kernel=kernel,
                                               scalekernel=scalekernel,
                                               edge=FALSE),
                                          list(...)))
               denominator <-
                 do.call(density.ppp,
                         resolve.defaults(list(x=X,
                                               at="pixels",
                                               weights = weights,
                                               sigma=sigma,
                                               varcov=varcov,
                                               kernel=kernel,
                                               scalekernel=scalekernel,
                                               edge=FALSE),
                                          list(...)))
               result <- eval.im(numerator/denominator)
               ## trap small values of denominator
               ## trap NaN and +/- Inf values of result, but not NA
               eps <- .Machine$double.eps
               nbg <- eval.im(is.infinite(result)
                              | is.nan(result)
                              | (denominator < eps))
               if(any(as.matrix(nbg), na.rm=TRUE)) {
                 warning(paste("Numerical underflow detected:",
                               "sigma is probably too small"))
                 ## l'Hopital's rule
                 distX <- distmap(X, xy=numerator)
                 whichnn <- attr(distX, "index")
                 nnvalues <- eval.im(values[whichnn])
                 result[nbg] <- nnvalues[nbg]
               }
               attr(result, "warnings") <- attr(numerator, "warnings")
             })
    }
  } else {
    ## ......... data frame of marks ..................
    ## detect constant columns
    ra <- apply(marx, 2, range, na.rm=TRUE)
    isconst <- (apply(ra, 2, diff) == 0)
    if(anyisconst <- any(isconst)) {
      oldmarx <- marx
#      oldX <- X
      marx <- marx[, !isconst]
      X <- X %mark% marx
    }
    if(any(!isconst)) {
      ## compute denominator
      denominator <-
        do.call(density.ppp,
                resolve.defaults(list(x=X,
                                      at=at,
                                      weights = weights,
                                      sigma=sigma, varcov=varcov,
                                      kernel=kernel,
                                      scalekernel=scalekernel,
                                      edge=FALSE),
                                 list(...)))
      ## compute numerator for each column of marks
      marx.weights <- if(weightsgiven) marx * weights else marx
      numerators <-
        do.call(density.ppp,
                resolve.defaults(list(x=X,
                                      at=at,
                                      weights = marx.weights,
                                      sigma=sigma, varcov=varcov,
                                      kernel=kernel,
                                      scalekernel=scalekernel,
                                      edge=FALSE),
                                 list(...)))
      uhoh <- attr(numerators, "warnings")
      ## calculate ratios
      switch(at,
             points={
               if(is.null(uhoh)) {
                 ## numerators is a matrix (or may have dropped to vector)
                 if(!is.matrix(numerators))
                   numerators <- matrix(numerators, ncol=1)
                 ratio <- numerators/denominator
                 if(any(badpoints <- matrowany(!is.finite(ratio)))) {
                   whichnnX <- nnwhich(X)
                   ratio[badpoints,] <-
                     as.matrix(marx[whichnnX[badpoints], , drop=FALSE])
                 }
               } else {
                 warning("returning original values")
                 ratio <- marx
               }
               result <- as.data.frame(ratio)
               colnames(result) <- colnames(marx)
             },
             pixels={
               ## numerators is a list of images (or may have dropped to 'im')
               if(is.im(numerators))
                 numerators <- list(numerators)
               ratio <- lapply(numerators, "/", e2=denominator)
               if(!is.null(uhoh)) {
                 ## compute nearest neighbour map on same raster
                 distX <- distmap(X, xy=denominator)
                 whichnnX <- attr(distX, "index")
                 ## fix images
                 for(j in 1:length(ratio)) {
                   ratj <- ratio[[j]]
                   valj <- marx[,j]
                   ratio[[j]] <-
                     eval.im(ifelseXY(is.finite(ratj), ratj, valj[whichnnX]))
                 }
                 attr(ratio, "warnings") <- uhoh
               }
               result <- as.solist(ratio)
               names(result) <- colnames(marx)
             })
    } else result <- NULL 
    if(anyisconst) {
      partresult <- result
      switch(at,
             points = {
               nX <- npoints(X)
               result <- matrix(, nX, ncol(oldmarx))
               if(length(partresult) > 0)
                 result[,!isconst] <- as.matrix(partresult)
               result[,isconst] <- rep(ra[1,isconst], each=nX)
               colnames(result) <- colnames(oldmarx)
             },
             pixels = {
               result <- vector(mode="list", length=ncol(oldmarx))
               if(length(partresult) > 0) {
                 result[!isconst] <- partresult
                 M <- as.owin(partresult[[1]])
               } else {
                 M <- do.call.matched(as.mask, list(w=as.owin(X), ...))
               }
               result[isconst] <- lapply(ra[1, isconst], as.im, W=M)
               result <- as.solist(result)
               names(result) <- colnames(oldmarx)
             })
    }
  }
  ## wrap up
  attr(result, "warnings") <-
    unlist(lapply(result, attr, which="warnings"))
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  return(result)
}


smoothpointsEngine <- function(x, values, sigma, ...,
                               kernel="gaussian", 
                               scalekernel=is.character(kernel),
                               weights=NULL, varcov=NULL,
                               leaveoneout=TRUE,
                               sorted=FALSE, cutoff=NULL) {
  debugging <- spatstat.options("developer")
  stopifnot(is.logical(leaveoneout))

  if(!is.null(dim(values)))
    stop("Internal error: smoothpointsEngine does not support multidimensional values")
  
  #' detect constant values
  if(diff(range(values, na.rm=TRUE)) == 0) { 
    result <- values
    attr(result, "sigma") <- sigma
    attr(result, "varcov") <- varcov
    return(result)
  }

  validate2Dkernel(kernel)
  if(is.character(kernel)) kernel <- match2DkernelName(kernel)
  isgauss <- identical(kernel, "gaussian")

  ## Handle weights that are meant to be null
  if(length(weights) == 0)
     weights <- NULL

  ## infinite bandwidth
  if(bandwidth.is.infinite(sigma)) {
    #' uniform estimate
    nX <- npoints(x)
    if(is.null(weights)) weights <- rep(1, nX)
    wtval <- weights * values
    totwt <- sum(weights)
    totwtval <- sum(wtval) 
    denominator <- rep(totwt, nX)
    numerator <- rep(totwtval, nX)
    if(leaveoneout) {
      numerator <- numerator - wtval
      denominator <- denominator - weights
    }
    result <- numerator/denominator
    return(result)
  }
  
  ## cutoff distance (beyond which the kernel value is treated as zero)
  ## NB: input argument 'cutoff' is either NULL or
  ##     an absolute distance (if scalekernel=FALSE)
  ##     a number of standard deviations (if scalekernel=TRUE)
  cutoff <- cutoff2Dkernel(kernel, sigma=sigma, varcov=varcov,
                           scalekernel=scalekernel, cutoff=cutoff,
                           fatal=TRUE)
  ## cutoff is now an absolute distance
  if(debugging)
    cat(paste("cutoff=", cutoff, "\n"))
  
  # detect very small bandwidth
  nnd <- nndist(x)
  nnrange <- range(nnd)
  if(cutoff < nnrange[1]) {
    if(leaveoneout && (npoints(x) > 1)) {
      warning("Very small bandwidth; values of nearest neighbours returned")
      result <- values[nnwhich(x)]
    } else {
      warning("Very small bandwidth; original values returned")
      result <- values
    }
    attr(result, "sigma") <- sigma
    attr(result, "varcov") <- varcov
    attr(result, "warnings") <- "underflow"
    return(result)
  }

  if(leaveoneout) {
    # ensure cutoff includes at least one point
    cutoff <- max(1.1 * nnrange[2], cutoff)
  }

  sd <- if(is.null(varcov)) sigma else sqrt(max(eigen(varcov)$values))
  
  if(isgauss &&
     spatstat.options("densityTransform") &&
     spatstat.options("densityC")) {
    ## .................. experimental C code .....................
    if(debugging)
      cat('Using experimental code!\n')
    npts <- npoints(x)
    result <- numeric(npts)
    ## transform to standard coordinates
    xx <- x$x
    yy <- x$y
    if(is.null(varcov)) {
      xx <- xx/(sqrt(2) * sigma)
      yy <- yy/(sqrt(2) * sigma)
    } else {
      Sinv <- solve(varcov)
      xy <- cbind(xx, yy) %*% matrixsqrt(Sinv/2)
      xx <- xy[,1]
      yy <- xy[,2]
      sorted <- FALSE
    }
    ## cutoff in standard coordinates
    cutoff <- cutoff/(sqrt(2) * sd)
    ## sort into increasing order of x coordinate (required by C code)
    if(!sorted) {
      oo <- fave.order(xx)
      xx <- xx[oo]
      yy <- yy[oo]
      vv <- values[oo]
    } else {
      vv <- values
    }
    if(is.null(weights)) {
      zz <- .C("Gsmoopt",
               nxy     = as.integer(npts),
               x       = as.double(xx),
               y       = as.double(yy),
               v       = as.double(vv),
               self    = as.integer(!leaveoneout),
               rmaxi   = as.double(cutoff),
               result  = as.double(double(npts)),
               PACKAGE = "spatstat")
      if(sorted) result <- zz$result else result[oo] <- zz$result
    } else {
      wtsort <- weights[oo]
      zz <- .C("Gwtsmoopt",
               nxy     = as.integer(npts),
               x       = as.double(xx),
               y       = as.double(yy),
               v       = as.double(vv),
               self    = as.integer(!leaveoneout),
               rmaxi   = as.double(cutoff),
               weight  = as.double(wtsort),
               result  = as.double(double(npts)),
               PACKAGE = "spatstat")
      if(sorted) result <- zz$result else result[oo] <- zz$result
    }
    if(any(nbg <- (is.infinite(result) | is.nan(result)))) {
      # NaN or +/-Inf can occur if bandwidth is small
      # Use mark of nearest neighbour (by l'Hopital's rule)
      result[nbg] <- values[nnwhich(x)[nbg]]
    }
  } else if(isgauss && spatstat.options("densityC")) {
    # .................. C code ...........................
    if(debugging)
      cat('Using standard code.\n')
    npts <- npoints(x)
    result <- numeric(npts)
    # sort into increasing order of x coordinate (required by C code)
    if(sorted) {
      xx <- x$x
      yy <- x$y
      vv <- values
    } else {
      oo <- fave.order(x$x)
      xx <- x$x[oo]
      yy <- x$y[oo]
      vv <- values[oo]
    }
    if(is.null(varcov)) {
      # isotropic kernel
      if(is.null(weights)) {
        zz <- .C("smoopt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 v       = as.double(vv),
                 self    = as.integer(!leaveoneout),
                 rmaxi   = as.double(cutoff),
                 sig     = as.double(sd),
                 result  = as.double(double(npts)),
                 PACKAGE = "spatstat")
        if(sorted) result <- zz$result else result[oo] <- zz$result
      } else {
        wtsort <- weights[oo]
        zz <- .C("wtsmoopt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 v       = as.double(vv),
                 self    = as.integer(!leaveoneout),
                 rmaxi   = as.double(cutoff),
                 sig     = as.double(sd),
                 weight  = as.double(wtsort),
                 result  = as.double(double(npts)),
                 PACKAGE = "spatstat")
        if(sorted) result <- zz$result else result[oo] <- zz$result
      }
    } else {
      # anisotropic kernel
      Sinv <- solve(varcov)
      flatSinv <- as.vector(t(Sinv))
      if(is.null(weights)) {
        zz <- .C("asmoopt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 v       = as.double(vv),
                 self    = as.integer(!leaveoneout),
                 rmaxi   = as.double(cutoff),
                 sinv    = as.double(flatSinv),
                 result  = as.double(double(npts)),
                 PACKAGE = "spatstat")
        if(sorted) result <- zz$result else result[oo] <- zz$result
      } else {
        wtsort <- weights[oo]
        zz <- .C("awtsmoopt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 v       = as.double(vv),
                 self    = as.integer(!leaveoneout),
                 rmaxi   = as.double(cutoff),
                 sinv    = as.double(flatSinv),
                 weight  = as.double(wtsort),
                 result  = as.double(double(npts)),
                 PACKAGE = "spatstat")
        if(sorted) result <- zz$result else result[oo] <- zz$result
      }
    }
    if(any(nbg <- (is.infinite(result) | is.nan(result)))) {
      # NaN or +/-Inf can occur if bandwidth is small
      # Use mark of nearest neighbour (by l'Hopital's rule)
      result[nbg] <- values[nnwhich(x)[nbg]]
    }
  } else {
    #' Either a non-Gaussian kernel or using older, partly interpreted code
    #' compute weighted densities
    if(is.null(weights)) {
      # weights are implicitly equal to 1
      numerator <- do.call(density.ppp,
                         resolve.defaults(list(x=x, at="points",
                                               weights = values,
                                               sigma=sigma,
                                               varcov=varcov,
                                               leaveoneout=leaveoneout,
                                               sorted=sorted,
                                               kernel=kernel,
                                               scalekernel=scalekernel),
                                          list(...),
                                          list(edge=FALSE)))
      denominator <- do.call(density.ppp,
                             resolve.defaults(list(x=x, at="points",
                                                   sigma=sigma,
                                                   varcov=varcov,
                                                   leaveoneout=leaveoneout,
                                                   sorted=sorted,
                                                   kernel=kernel,
                                                   scalekernel=scalekernel),
                                              list(...),
                                              list(edge=FALSE)))
    } else {
      numerator <- do.call(density.ppp,
                           resolve.defaults(list(x=x, at="points",
                                                 weights = values * weights,
                                                 sigma=sigma,
                                                 varcov=varcov,
                                                 leaveoneout=leaveoneout,
                                                 sorted=sorted,
                                                 kernel=kernel,
                                                 scalekernel=scalekernel),
                                            list(...),
                                            list(edge=FALSE)))
      denominator <- do.call(density.ppp,
                             resolve.defaults(list(x=x, at="points",
                                                   weights = weights,
                                                   sigma=sigma,
                                                   varcov=varcov,
                                                   leaveoneout=leaveoneout,
                                                   sorted=sorted,
                                                   kernel=kernel,
                                                   scalekernel=scalekernel),
                                              list(...),
                                              list(edge=FALSE)))
    }
    if(is.null(uhoh <- attr(numerator, "warnings"))) {
      result <- numerator/denominator
      result <- ifelseXB(is.finite(result), result, NA_real_)
    } else {
      warning("returning original values")
      result <- values
      attr(result, "warnings") <- uhoh
    }
  }
  # pack up and return
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  return(result)
}


markmean <- function(X, ...) {
  stopifnot(is.marked(X))
  Y <- Smooth(X, ...)
  return(Y)
}

markvar  <- function(X, sigma=NULL, ..., weights=NULL, varcov=NULL) {
  stopifnot(is.marked(X))
  if(is.expression(weights)) 
    weights <- eval(weights, envir=as.data.frame(X), enclos=parent.frame())
  E1 <- Smooth(X, sigma=sigma, varcov=varcov, weights=weights, ...)
  X2 <- X %mark% marks(X)^2
  ## ensure smoothing bandwidth is the same!
  sigma <- attr(E1, "sigma")
  varcov <- attr(E1, "varcov")
  E2 <- Smooth(X2, sigma=sigma, varcov=varcov, weights=weights, ...)
  V <- eval.im(E2 - E1^2)
  return(V)
}

bw.smoothppp <- function(X, nh=spatstat.options("n.bandwidth"),
                         hmin=NULL, hmax=NULL, warn=TRUE,
                         kernel="gaussian") {
  stopifnot(is.ppp(X))
  stopifnot(is.marked(X))
  X <- coerce.marks.numeric(X)
  # rearrange in ascending order of x-coordinate (for C code)
  X <- X[fave.order(X$x)]
  #
  marx <- marks(X)
  dimmarx <- dim(marx)
  if(!is.null(dimmarx))
    marx <- as.matrix(as.data.frame(marx))
  # determine a range of bandwidth values
#  n <- npoints(X)
  if(is.null(hmin) || is.null(hmax)) {
    W <- Window(X)
#    a <- area(W)
    d <- diameter(as.rectangle(W))
    # Stoyan's rule of thumb 
    stoyan <- bw.stoyan(X)
    # rule of thumb based on nearest-neighbour distances
    nnd <- nndist(X)
    nnd <- nnd[nnd > 0]
    if(is.null(hmin)) {
      hmin <- max(1.1 * min(nnd), stoyan/5)
      hmin <- min(d/8, hmin)
    }
    if(is.null(hmax)) {
      hmax <- max(stoyan * 20, 3 * mean(nnd), hmin * 2)
      hmax <- min(d/2, hmax)
    }
  } else stopifnot(hmin < hmax)
  #
  h <- geomseq(from=hmin, to=hmax, length.out=nh)
  cv <- numeric(nh)
  # 
  # compute cross-validation criterion
  for(i in seq_len(nh)) {
    yhat <- Smooth(X, sigma=h[i], at="points", leaveoneout=TRUE,
                   kernel=kernel, 
                   sorted=TRUE)
    if(!is.null(dimmarx))
      yhat <- as.matrix(as.data.frame(yhat))
    cv[i] <- mean((marx - yhat)^2)
  }

  # optimize
  iopt <- which.min(cv)
#  hopt <- h[iopt]
  #
  if(warn && (iopt == nh || iopt == 1)) 
    warning(paste("Cross-validation criterion was minimised at",
                  if(iopt == 1) "left-hand" else "right-hand",
                  "end of interval",
                  paste(prange(signif(c(hmin, hmax), 3)), ";", sep=""),
                  "use arguments hmin, hmax to specify a wider interval"),
            call.=FALSE)
  #
  result <- bw.optim(cv, h, iopt,
                     hname="sigma",
                     creator="bw.smoothppp",
                     criterion="Least Squares Cross-Validation",
                     unitname=unitname(X))
  return(result)
}

smoothcrossEngine <- function(Xdata, Xquery, values, sigma, ...,
                              weights=NULL, varcov=NULL,
                              kernel="gaussian", 
                              scalekernel=is.character(kernel),
                              sorted=FALSE,
                              cutoff=NULL) {

  validate2Dkernel(kernel)
  if(is.character(kernel)) kernel <- match2DkernelName(kernel)
  isgauss <- identical(kernel, "gaussian") && scalekernel

  if(!is.null(dim(weights)))
    stop("weights must be a vector")

  ndata <- npoints(Xdata)
  nquery <- npoints(Xquery)
  
  if(nquery == 0 || ndata == 0) {
    if(is.null(dim(values))) return(rep(NA_real_, nquery))
    nuttin <- matrix(NA_real_, nrow=nquery, ncol=ncol(values))
    colnames(nuttin) <- colnames(values)
    return(nuttin)
  }

  # validate weights
  if(is.matrix(values) || is.data.frame(values)) {
    k <- ncol(values)
    stopifnot(nrow(values) == npoints(Xdata))
    values <- as.data.frame(values)
  } else {
    k <- 1L
    stopifnot(length(values) == npoints(Xdata) || length(values) == 1)
    if(length(values) == 1L)
      values <- rep(values, ndata)
  }

  ## infinite bandwidth
  if(bandwidth.is.infinite(sigma)) {
    #' uniform estimate
    if(is.null(weights)) weights <- rep(1, ndata)
    single <- is.null(dim(values))
    wtval <- weights * values 
    totwt <- sum(weights)
    totwtval <- if(single) sum(wtval) else colSums(wtval)
    denominator <- rep(totwt, nquery)
    numerator <- rep(totwtval, each=nquery)
    if(!single) numerator <- matrix(numerator, nrow=nquery)
    result <- numerator/denominator
    if(!single)
      colnames(result) <- colnames(values)
    return(result)
  }
  
  ## cutoff distance (beyond which the kernel value is treated as zero)
  ## NB: input argument 'cutoff' is either NULL or
  ##     an absolute distance (if scalekernel=FALSE)
  ##     a number of standard deviations (if scalekernel=TRUE)
  cutoff.orig <- cutoff
  cutoff <- cutoff2Dkernel(kernel, sigma=sigma, varcov=varcov,
                           scalekernel=scalekernel, cutoff=cutoff,
                           fatal=TRUE)
  ## cutoff is now an absolute distance

  ## detect very small bandwidth
  nnc <- nncross(Xquery, Xdata)
  if(cutoff < min(nnc$dist)) {
    if(ndata > 1) {
      warning("Very small bandwidth; values of nearest neighbours returned")
      nw <- nnc$which
      result <- if(k == 1) values[nw] else values[nw,,drop=FALSE]
    } else {
      warning("Very small bandwidth; original values returned")
      result <- values
    }
    attr(result, "sigma") <- sigma
    attr(result, "varcov") <- varcov
    attr(result, "warnings") <- "underflow"
    return(result)
  }
  
  ## Handle weights that are meant to be null
  if(length(weights) == 0)
     weights <- NULL

  if(!isgauss) {
    ## .................. non-Gaussian kernel ........................
    close <- crosspairs(Xdata, Xquery, cutoff)
    kerij <- evaluate2Dkernel(kernel, close$dx, close$dy,
                            sigma=sigma, varcov=varcov,
                            scalekernel=scalekernel, ...)
    ## sum the (weighted) contributions
    i <- close$i # data point
    j <- close$j # query point
    jfac <- factor(j, levels=seq_len(nquery))
    wkerij <- if(is.null(weights)) kerij else kerij * weights[i]
    denominator <- tapplysum(wkerij, list(jfac))
    if(k == 1L) {
      contribij <- wkerij * values[i]
      numerator <- tapplysum(contribij, list(jfac))
      result <- numerator/denominator
    } else {
      for(kk in 1:k) {
        contribij <- wkerij * values[i, kk]
        numeratorkk <- tapplysum(contribij, list(jfac))
        result[,kk] <- numeratorkk/denominator
      }
    }
    ## trap bad values
    if(any(nbg <- (is.infinite(result) | is.nan(result)))) {
      ## NaN or +/-Inf can occur if bandwidth is small
      ## Use value at nearest neighbour (by l'Hopital's rule)
      nnw <- nnc$which
      if(k == 1L) {
        result[nbg] <- values[nnw[nbg]]
      } else {
        bad <- which(nbg, arr.ind=TRUE)
        badrow <- bad[,"row"]
        badcol <- bad[,"col"]
        result[nbg] <- values[cbind(nnw[badrow], badcol)]
      }
    }
    attr(result, "sigma") <- sigma
    attr(result, "varcov") <- varcov
    return(result)
  } 

  ## .................. Gaussian kernel henceforth ........................
  
  ## handle multiple columns of values
  if(is.matrix(values) || is.data.frame(values)) {
    k <- ncol(values)
    stopifnot(nrow(values) == npoints(Xdata))
    values <- as.data.frame(values)
    result <- matrix(, nquery, k)
    colnames(result) <- colnames(values)
    if(!sorted) {
      ood <- fave.order(Xdata$x)
      Xdata <- Xdata[ood]
      values <- values[ood, ]
      ooq <- fave.order(Xquery$x)
      Xquery <- Xquery[ooq]
    }
    for(j in 1:k) 
      result[,j] <- smoothcrossEngine(Xdata, Xquery, values[,j],
                                      sigma=sigma, varcov=varcov,
                                      weights=weights,
                                      kernel=kernel, scalekernel=scalekernel,
                                      cutoff=cutoff.orig,
                                      sorted=TRUE, ...)
    if(!sorted) {
      sortresult <- result
      result[ooq,] <- sortresult
    }
    attr(result, "sigma") <- sigma
    attr(result, "varcov") <- varcov
    return(result)
  }

  ## values must be a vector
  stopifnot(length(values) == npoints(Xdata) || length(values) == 1)
  if(length(values) == 1) values <- rep(values, ndata)

  result <- numeric(nquery) 
  ## coordinates and values
  xq <- Xquery$x
  yq <- Xquery$y
  xd <- Xdata$x
  yd <- Xdata$y
  vd <- values
  if(!sorted) {
    ## sort into increasing order of x coordinate (required by C code)
    ooq <- fave.order(Xquery$x)
    xq <- xq[ooq]
    yq <- yq[ooq]
    ood <- fave.order(Xdata$x)
    xd <- xd[ood]
    yd <- yd[ood]
    vd <- vd[ood] 
  }
  sd <- if(is.null(varcov)) sigma else sqrt(min(eigen(varcov)$values))
  if(is.null(varcov)) {
    ## isotropic kernel
    if(is.null(weights)) {
      zz <- .C("crsmoopt",
               nquery      = as.integer(nquery),
               xq      = as.double(xq),
               yq      = as.double(yq),
               ndata   = as.integer(ndata),
               xd      = as.double(xd),
               yd      = as.double(yd),
               vd      = as.double(vd),
               rmaxi   = as.double(cutoff),
               sig     = as.double(sd),
               result  = as.double(double(nquery)),
               PACKAGE = "spatstat")
      if(sorted) result <- zz$result else result[ooq] <- zz$result
    } else {
      wtsort <- if(sorted) weights else weights[ood]
      zz <- .C("wtcrsmoopt",
               nquery      = as.integer(nquery),
               xq      = as.double(xq),
               yq      = as.double(yq),
               ndata   = as.integer(ndata),
               xd      = as.double(xd),
               yd      = as.double(yd),
               vd      = as.double(vd),
               wd      = as.double(wtsort),
               rmaxi   = as.double(cutoff),
               sig     = as.double(sd),
               result  = as.double(double(nquery)),
               PACKAGE = "spatstat")
        if(sorted) result <- zz$result else result[ooq] <- zz$result
      }
    } else {
      # anisotropic kernel
      Sinv <- solve(varcov)
      flatSinv <- as.vector(t(Sinv))
      if(is.null(weights)) {
        zz <- .C("acrsmoopt",
                 nquery      = as.integer(nquery),
                 xq      = as.double(xq),
                 yq      = as.double(yq),
                 ndata   = as.integer(ndata),
                 xd      = as.double(xd),
                 yd      = as.double(yd),
                 vd      = as.double(vd),
                 rmaxi   = as.double(cutoff),
                 sinv    = as.double(flatSinv),
                 result  = as.double(double(nquery)),
                 PACKAGE = "spatstat")
        if(sorted) result <- zz$result else result[ooq] <- zz$result
      } else {
        wtsort <- if(sorted) weights else weights[ood]
        zz <- .C("awtcrsmoopt",
                 nquery      = as.integer(nquery),
                 xq      = as.double(xq),
                 yq      = as.double(yq),
                 ndata   = as.integer(ndata),
                 xd      = as.double(xd),
                 yd      = as.double(yd),
                 vd      = as.double(vd),
                 wd      = as.double(wtsort),
                 rmaxi   = as.double(cutoff),
                 sinv    = as.double(flatSinv),
                 result  = as.double(double(nquery)),
                 PACKAGE = "spatstat")
        if(sorted) result <- zz$result else result[ooq] <- zz$result
      }
    }
    if(any(nbg <- (is.infinite(result) | is.nan(result)))) {
      # NaN or +/-Inf can occur if bandwidth is small
      # Use mark of nearest neighbour (by l'Hopital's rule)
      result[nbg] <- values[nnc$which[nbg]]
    }
  # pack up and return
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  return(result)
}

ExpSmoothLog <- function(X, ..., at=c("pixels", "points"), weights=NULL) {
  verifyclass(X, "ppp")
  at <- match.arg(at)
  if(!is.null(weights)) 
    check.nvector(weights, npoints(X))
  X <- coerce.marks.numeric(X)
  marx <- marks(X)
  d <- dim(marx)
  if(!is.null(d) && d[2] > 1) {
    switch(at,
           points = {
             Z <- lapply(unstack(X), ExpSmoothLog, ...,
                         at=at, weights=weights)
             Z <- do.call(data.frame, Z)
           },
           pixels = {
             Z <- solapply(unstack(X), ExpSmoothLog, ...,
                           at=at, weights=weights)
           })
    return(Z)
  }
  # vector or single column of numeric marks
  v <- as.numeric(marx)
  vmin <- min(v)
  if(vmin < 0) stop("Negative values in geometric mean smoothing",
                       call.=FALSE)
  Y <- X %mark% log(v)
  if(vmin > 0) {
    Z <- Smooth(Y, ..., at=at, weights=weights)
  } else {
    yok <- is.finite(marks(Y))
    YOK <- Y[yok]
    weightsOK <- if(is.null(weights)) NULL else weights[yok]
    switch(at,
           points = {
             Z <- rep(-Inf, npoints(X))
             Z[yok] <- Smooth(YOK, ..., at=at, weights=weightsOK)
           },
           pixels = {
             isfinite <- nnmark(Y %mark% yok, ...)
             support <- solutionset(isfinite)
             Window(YOK) <- support
             Z <- as.im(-Inf, W=Window(Y), ...)
             Z[support] <- Smooth(YOK, ..., at=at, weights=weightsOK)[]
           })
  }
  return(exp(Z))
}
