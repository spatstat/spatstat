#
#  smooth.ppp.R
#
#  Smooth the marks of a point pattern
# 
#  $Revision: 1.4 $  $Date: 2013/05/01 08:03:17 $
#

smooth.ppp <- function(X, ..., weights=rep(1, npoints(X)), at="pixels") {
  verifyclass(X, "ppp")
  if(!is.marked(X, dfok=TRUE))
    stop("X should be a marked point pattern")
  at <- pickoption("output location type", at,
                   c(pixels="pixels",
                     points="points"))
  # determine smoothing parameters
  ker <- resolve.2D.kernel(..., x=X, bwfun=bw.smoothppp)
  sigma <- ker$sigma
  varcov <- ker$varcov
  #
  if(weightsgiven <- !missing(weights)) {
    check.nvector(weights, npoints(X)) 
    # rescale weights to avoid numerical gremlins
    weights <- weights/median(abs(weights))
  }
  # get marks
  marx <- marks(X, dfok=TRUE)
  #
  if(!is.data.frame(marx)) {
    # ........ vector of marks ...................
    values <- marx
    if(is.factor(values)) {
      warning("Factor valued marks were converted to integers")
      values <- as.numeric(values)
    }
    switch(at,
           points={
             if(!weightsgiven)
               weights <- NULL
             result <-
               do.call("smoothpointsEngine",
                       resolve.defaults(list(x=X),
                                        list(values=values, weights=weights),
                                        list(sigma=sigma, varcov=varcov),
                                        list(...)))
           },
           pixels={
             numerator <-
               do.call("density.ppp",
                       resolve.defaults(list(x=X, at="pixels"),
                                        list(weights = values * weights),
                                        list(sigma=sigma, varcov=varcov),
                                        list(...),
                                        list(edge=FALSE)))
             denominator <-
               do.call("density.ppp",
                       resolve.defaults(list(x=X, at="pixels"),
                                        list(weights = weights),
                                        list(sigma=sigma, varcov=varcov),
                                        list(...),
                                        list(edge=FALSE)))
             result <- eval.im(numerator/denominator)
             # trap small values of denominator
             # trap NaN and +/- Inf values of result, but not NA
             eps <- .Machine$double.eps
             nbg <- eval.im(is.infinite(result)
                            | is.nan(result)
                            | (denominator < eps))
             if(any(as.matrix(nbg), na.rm=TRUE)) {
               warning("Numerical underflow detected: sigma is probably too small")
               # l'Hopital's rule
               distX <- distmap(X, xy=numerator)
               whichnn <- attr(distX, "index")
               nnvalues <- eval.im(values[whichnn])
               result[nbg] <- nnvalues[nbg]
             }
             attr(result, "warnings") <- attr(numerator, "warnings")
           })
  } else {
    # ......... data frame of marks ..................
    nmarx <- ncol(marx)
    result <- list()
    # compute denominator
    denominator <-
      do.call("density.ppp",
              resolve.defaults(list(x=X, at=at),
                               list(weights = weights),
                               list(sigma=sigma, varcov=varcov),
                               list(...),
                               list(edge=FALSE)))
    switch(at,
           pixels = {
             # compute nearest neighbour map on same raster
             distX <- distmap(X, xy=denominator)
             whichnnX <- attr(distX, "index")
           },
           points = {
             whichnnX <- nnwhich(X)
           })
    # smooth each column of marks in turn
    for(j in 1:nmarx) {
      values <- marx[,j] 
      if(is.factor(values)) {
        warning(paste("Factor valued marks in column", j,
                      "were converted to integers"))
        values <- as.numeric(values)
      }
      # compute j-th numerator
      numerator <-
        do.call("density.ppp",
                resolve.defaults(list(x=X, at=at),
                                 list(weights = values * weights),
                                 list(sigma=sigma, varcov=varcov),
                                 list(...),
                                 list(edge=FALSE)))
      switch(at,
             points={
               if(is.null(uhoh <- attr(numerator, "warnings"))) {
                 ratio <- numerator/denominator
                 ratio <- ifelseXY(is.finite(ratio), ratio, values[whichnnX])
               } else {
                 warning("returning original values")
                 ratio <- values
                 attr(ratio, "warnings") <- uhoh
               }
             },
             pixels={
               ratio <- eval.im(numerator/denominator)
               ratio <- eval.im(ifelseXY(is.finite(ratio), ratio,
                                         values[whichnnX]))
               attr(ratio, "warnings") <- attr(numerator, "warnings")
           })
    # store results
      result[[j]] <- ratio
    }
    # wrap up
    names(result) <- colnames(marx)
    result <- switch(at,
                     pixels=as.listof(result),
                     points=as.data.frame(result))
    attr(result, "warnings") <-
      unlist(lapply(result, function(x){ attr(x, "warnings") }))
  }
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  return(result)
}

smoothpointsEngine <- function(x, values, sigma, ...,
                               weights=NULL, varcov=NULL,
                               leaveoneout=TRUE,
                               sorted=FALSE) {
  stopifnot(is.logical(leaveoneout))
  if(is.null(varcov)) {
    const <- 1/(2 * pi * sigma^2)
  } else {
    detSigma <- det(varcov)
    Sinv <- solve(varcov)
    const <- 1/(2 * pi * sqrt(detSigma))
  }
  # Contributions from pairs of distinct points
  # closer than 8 standard deviations
  sd <- if(is.null(varcov)) sigma else sqrt(sum(diag(varcov)))
  cutoff <- 8 * sd
  
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
                 result   = as.double(double(npts)),
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
                 result   = as.double(double(npts)),
                 PACKAGE = "spatstat")
        if(sorted) result <- zz$result else result[oo] <- zz$result
      }
    } else {
      # anisotropic kernel
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
                 result   = as.double(double(npts)),
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
                 result   = as.double(double(npts)),
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


markmean <- function(X, ...) { smooth.ppp(X, ...) }

markvar  <- function(X, ...) {
  if(!is.marked(X, dfok=FALSE))
    stop("X should have (one column of) marks")
  E1 <- smooth.ppp(X, ...)
  X2 <- X %mark% marks(X)^2
  E2 <- smooth.ppp(X2, ...)
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
  n <- npoints(X)
  if(is.null(hmin) || is.null(hmax)) {
    W <- as.owin(X)
    a <- area.owin(W)
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
    yhat <- smooth.ppp(X, sigma=h[i], at="points", leaveoneout=TRUE,
                       sorted=TRUE)
    cv[i] <- mean((marx - yhat)^2)
  }

  # optimize
  iopt <- which.min(cv)
  hopt <- h[iopt]
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
