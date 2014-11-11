#
#  smooth.ppp.R
#
#  Smooth the marks of a point pattern
# 
#  $Revision: 1.22 $  $Date: 2014/10/24 00:22:30 $
#

smooth.ppp <- function(X, ..., weights=rep(1, npoints(X)), at="pixels") {
  .Deprecated("Smooth.ppp", package="spatstat",
    msg="smooth.ppp is deprecated: use the generic Smooth with a capital S")
  Smooth(X, ..., weights=weights, at=at)
}

Smooth <- function(X, ...) {
  UseMethod("Smooth")
}

Smooth.ppp <- function(X, sigma=NULL, ...,
                       weights=rep(1, npoints(X)), at="pixels",
                       edge=TRUE, diggle=FALSE) {
  verifyclass(X, "ppp")
  if(!is.marked(X, dfok=TRUE))
    stop("X should be a marked point pattern")
  X <- coerce.marks.numeric(X)
  at <- pickoption("output location type", at,
                   c(pixels="pixels",
                     points="points"))
  ## determine smoothing parameters
  ker <- resolve.2D.kernel(sigma=sigma, ...,
                           x=X, bwfun=bw.smoothppp, allow.zero=TRUE)
  sigma <- ker$sigma
  varcov <- ker$varcov
  ## Diggle's edge correction?
  if(diggle && !edge) warning("Option diggle=TRUE overridden by edge=FALSE")
  diggle <- diggle && edge
  ## 
  if(ker$cutoff < minnndist(X)) {
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
  ## weights
  weightsgiven <- !missing(weights) && !is.null(weights) && (length(weights)>0)
  if(weightsgiven) {
    check.nvector(weights, npoints(X))
  } else weights <- NULL
  if(diggle) {
    ## absorb Diggle edge correction into weights vector
    edg <- second.moment.calc(X, sigma, what="edge", ..., varcov=varcov)
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
                 do.call("smoothpointsEngine",
                         resolve.defaults(list(x=X,
                                               values=values, weights=weights,
                                               sigma=sigma, varcov=varcov,
                                               edge=FALSE),
                                          list(...)))
             },
             pixels={
               values.weights <- if(weightsgiven) values * weights else values
               numerator <-
                 do.call("density.ppp",
                         resolve.defaults(list(x=X,
                                               at="pixels",
                                               weights = values.weights,
                                               sigma=sigma, varcov=varcov,
                                               edge=FALSE),
                                          list(...)))
               denominator <-
                 do.call("density.ppp",
                         resolve.defaults(list(x=X,
                                               at="pixels",
                                               weights = weights,
                                               sigma=sigma,
                                               varcov=varcov,
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
        do.call("density.ppp",
                resolve.defaults(list(x=X,
                                      at=at,
                                      weights = weights,
                                      sigma=sigma, varcov=varcov,
                                      edge=FALSE),
                                 list(...)))
      ## compute numerator for each column of marks
      marx.weights <- if(weightsgiven) marx * weights else marx
      numerators <-
        do.call("density.ppp",
                resolve.defaults(list(x=X,
                                      at=at,
                                      weights = marx.weights,
                                      sigma=sigma, varcov=varcov,
                                      edge=FALSE),
                                 list(...)))
      uhoh <- attr(numerators, "warnings")
      ## calculate ratios
      switch(at,
             points={
               if(is.null(uhoh)) {
                 ## numerators is a matrix
                 ratio <- numerators/denominator
                 if(any(badpoints <- apply(!is.finite(ratio), 1, any))) {
                   whichnnX <- nnwhich(X)
                   ratio[badpoints,] <- marx[whichnnX[badpoints], ]
                 }
               } else {
                 warning("returning original values")
                 ratio <- marx
               }
               result <- as.data.frame(ratio)
               colnames(result) <- colnames(marx)
             },
             pixels={
               ratio <- lapply(numerators,
                               function(a,b) eval.im(a/b),
                               b=denominator)
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
               result <- as.listof(ratio)
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
                 result[,!isconst] <- partresult
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
               result <- as.listof(result)
               names(result) <- colnames(oldmarx)
             })
    }
  }
  ## wrap up
  attr(result, "warnings") <-
    unlist(lapply(result, function(x){ attr(x, "warnings") }))
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  return(result)
}


smoothpointsEngine <- function(x, values, sigma, ...,
                               weights=NULL, varcov=NULL,
                               leaveoneout=TRUE,
                               sorted=FALSE) {
  stopifnot(is.logical(leaveoneout))
#  if(is.null(varcov)) {
#    const <- 1/(2 * pi * sigma^2)
#  } else {
#    detSigma <- det(varcov)
#    Sinv <- solve(varcov)
#    const <- 1/(2 * pi * sqrt(detSigma))
#  }
  # detect constant values
  if(diff(range(values, na.rm=TRUE)) == 0) { 
    result <- values
    attr(result, "sigma") <- sigma
    attr(result, "varcov") <- varcov
    return(result)
  }
  # Contributions from pairs of distinct points
  # closer than 8 standard deviations
  sd <- if(is.null(varcov)) sigma else sqrt(sum(diag(varcov)))
  cutoff <- 8 * sd

  ## Handle weights that are meant to be null
  if(length(weights) == 0 || (!is.null(dim(weights)) && nrow(weights) == 0))
     weights <- NULL
     
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
  if(spatstat.options("densityC")) {
    # .................. new C code ...........................
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
                 result  = as.double(double(npts)))
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
                 result  = as.double(double(npts)))
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
                 result  = as.double(double(npts)))
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
                 result  = as.double(double(npts)))
        if(sorted) result <- zz$result else result[oo] <- zz$result
      }
    }
    if(any(nbg <- (is.infinite(result) | is.nan(result)))) {
      # NaN or +/-Inf can occur if bandwidth is small
      # Use mark of nearest neighbour (by l'Hopital's rule)
      result[nbg] <- values[nnwhich(x)[nbg]]
    }
  } else {
    # previous, partly interpreted code
    # compute weighted densities
    if(is.null(weights)) {
      # weights are implicitly equal to 1
      numerator <- do.call("density.ppp",
                         resolve.defaults(list(x=x, at="points"),
                                          list(weights = values),
                                          list(sigma=sigma, varcov=varcov),
                                          list(leaveoneout=leaveoneout),
                                          list(sorted=sorted),
                                          list(...),
                                          list(edge=FALSE)))
      denominator <- do.call("density.ppp",
                             resolve.defaults(list(x=x, at="points"),
                                              list(sigma=sigma, varcov=varcov),
                                              list(leaveoneout=leaveoneout),
                                              list(sorted=sorted),
                                              list(...),
                                              list(edge=FALSE)))
    } else {
      numerator <- do.call("density.ppp",
                           resolve.defaults(list(x=x, at="points"),
                                            list(weights = values * weights),
                                            list(sigma=sigma, varcov=varcov),
                                            list(leaveoneout=leaveoneout),
                                            list(sorted=sorted),
                                            list(...),
                                            list(edge=FALSE)))
      denominator <- do.call("density.ppp",
                             resolve.defaults(list(x=x, at="points"),
                                              list(weights = weights),
                                              list(sigma=sigma, varcov=varcov),
                                              list(leaveoneout=leaveoneout),
                                              list(sorted=sorted),
                                              list(...),
                                              list(edge=FALSE)))
    }
    if(is.null(uhoh <- attr(numerator, "warnings"))) {
      result <- numerator/denominator
      result <- ifelseXB(is.finite(result), result, NA)
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
                       hmin=NULL, hmax=NULL, warn=TRUE) {
  stopifnot(is.ppp(X))
  stopifnot(is.marked(X))
  # rearrange in ascending order of x-coordinate (for C code)
  X <- X[fave.order(X$x)]
  #
  marx <- marks(X)
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
  h <- exp(seq(from=log(hmin), to=log(hmax), length.out=nh))
  cv <- numeric(nh)
  # 
  # compute cross-validation criterion
  for(i in seq_len(nh)) {
    yhat <- Smooth(X, sigma=h[i], at="points", leaveoneout=TRUE,
                       sorted=TRUE)
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
                     criterion="Least Squares Cross-Validation")
  return(result)
}

smoothcrossEngine <- function(Xdata, Xquery, values, sigma, ...,
                              weights=NULL, varcov=NULL,
                              sorted=FALSE) {
#  if(is.null(varcov)) {
#    const <- 1/(2 * pi * sigma^2)
#  } else {
#    detSigma <- det(varcov)
#    Sinv <- solve(varcov)
#    const <- 1/(2 * pi * sqrt(detSigma))
#  }
  if(!is.null(dim(weights)))
    stop("weights must be a vector")

  if(npoints(Xquery) == 0 || npoints(Xdata) == 0) {
    if(is.null(dim(values))) return(rep(NA, npoints(Xquery)))
    nuttin <- matrix(NA, nrow=npoints(Xquery), ncol=ncol(values))
    colnames(nuttin) <- colnames(values)
    return(nuttin)
  }
  
  ## Contributions from pairs of distinct points
  ## closer than 8 standard deviations
  sd <- if(is.null(varcov)) sigma else sqrt(sum(diag(varcov)))
  cutoff <- 8 * sd

  ## detect very small bandwidth
  nnc <- nncross(Xquery, Xdata)
  if(cutoff < min(nnc$dist)) {
    if(npoints(Xdata) > 1) {
      warning("Very small bandwidth; values of nearest neighbours returned")
      nw <- nnc$which
      result <- if(is.null(dim(values))) values[nw] else values[nw,,drop=FALSE]
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
     
  ## handle multiple columns of values
  if(is.matrix(values) || is.data.frame(values)) {
    k <- ncol(values)
    stopifnot(nrow(values) == npoints(Xdata))
    values <- as.data.frame(values)
    result <- matrix(, npoints(Xdata), k)
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
                                      weights=weights, sorted=TRUE,
                                      ...)
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
  if(length(values) == 1) values <- rep(values, npoints(Xdata))

  ndata <- npoints(Xdata)
  nquery <- npoints(Xquery)
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
               result  = as.double(double(nquery)))
      if(sorted) result <- zz$result else result[ooq] <- zz$result
    } else {
      wtsort <- weights[ood]
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
               result  = as.double(double(nquery)))
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
                 result  = as.double(double(nquery)))
        if(sorted) result <- zz$result else result[ooq] <- zz$result
      } else {
        wtsort <- weights[ood]
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
                 result  = as.double(double(nquery)))
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

