#
#  density.ppp.R
#
#  Method for 'density' for point patterns
#
#  $Revision: 1.85 $    $Date: 2016/11/13 01:54:32 $
#

ksmooth.ppp <- function(x, sigma, ..., edge=TRUE) {
  .Deprecated("density.ppp", package="spatstat")
  density.ppp(x, sigma, ..., edge=edge)
}

density.ppp <- local({
  
density.ppp <- function(x, sigma=NULL, ...,
                        weights=NULL, edge=TRUE, varcov=NULL,
                        at="pixels", leaveoneout=TRUE,
                        adjust=1, diggle=FALSE, se=FALSE, 
                        kernel="gaussian",
                        scalekernel=is.character(kernel),
                        positive=FALSE) {
  verifyclass(x, "ppp")

  output <- pickoption("output location type", at,
                       c(pixels="pixels",
                         points="points"))

  if(!identical(kernel, "gaussian")) {
    validate2Dkernel(kernel)
    ## kernel is only partly implemented!
    if(se)
      stop("Standard errors are not implemented for non-Gaussian kernel")
    if(is.function(sigma) || (is.null(sigma) && is.null(varcov)))
      warning("Bandwidth selection will be based on Gaussian kernel")
  }
  
  ker <- resolve.2D.kernel(..., sigma=sigma, varcov=varcov, x=x, adjust=adjust)
  sigma <- ker$sigma
  varcov <- ker$varcov

  if(is.im(weights)) {
    weights <- safelookup(weights, x) # includes warning if NA
  } else if(is.expression(weights)) 
    weights <- eval(weights, envir=as.data.frame(x), enclos=parent.frame())
  if(length(weights) == 0 || (!is.null(dim(weights)) && nrow(weights) == 0))
    weights <- NULL 

  if(se) {
    # compute standard error
    SE <- denspppSEcalc(x, sigma=sigma, varcov=varcov,
                        ...,
                        weights=weights, edge=edge, at=output,
                        leaveoneout=leaveoneout, adjust=adjust,
                        diggle=diggle)
    if(positive) SE <- posify(SE)
  }
                         
  if(output == "points") {
    # VALUES AT DATA POINTS ONLY
    result <- densitypointsEngine(x, sigma,
                                  varcov=varcov,
                                  kernel=kernel,
                                  scalekernel=scalekernel,
                                  weights=weights, edge=edge,
                                  leaveoneout=leaveoneout,
                                  diggle=diggle, ...)
    if(!is.null(uhoh <- attr(result, "warnings"))) {
      switch(uhoh,
             underflow=warning("underflow due to very small bandwidth"),
             warning(uhoh))
    }
    ## constrain values to be positive
    if(positive) 
      result <- posify(result)
    if(se) 
      result <- list(estimate=result, SE=SE)
    return(result)
  }
  
  # VALUES AT PIXELS
  if(!edge) {
    # no edge correction
    edg <- NULL
    raw <- second.moment.calc(x, sigma, what="smooth", ...,
                              kernel=kernel,
                              scalekernel=scalekernel,
                              weights=weights, varcov=varcov)
    raw <- divide.by.pixelarea(raw) 
    smo <- raw
  } else if(!diggle) {
    # edge correction e(u)
    both <- second.moment.calc(x, sigma, what="smoothedge", ...,
                               kernel=kernel,
                               scalekernel=scalekernel,
                               weights=weights, varcov=varcov)
    raw <- divide.by.pixelarea(both$smooth)
    edg <- both$edge
    smo <- if(is.im(raw)) eval.im(raw/edg) else
           lapply(raw, divideimage, denom=edg)
  } else {
    # edge correction e(x_i)
    edg <- second.moment.calc(x, sigma, what="edge", ...,
                              scalekernel=scalekernel,
                              kernel=kernel, varcov=varcov)
    wi <- 1/safelookup(edg, x, warn=FALSE)
    wi[!is.finite(wi)] <- 0
    # edge correction becomes weight attached to points
    if(is.null(weights)) {
      newweights <- wi
    } else if(is.matrix(weights) || is.data.frame(weights)) {
      stopifnot(nrow(weights) == npoints(x))
      newweights <- weights * wi
    } else {
      stopifnot(length(weights) == npoints(x))
      newweights <- weights * wi
    }
    raw <- second.moment.calc(x, sigma, what="smooth", ...,
                              kernel=kernel,
                              scalekernel=scalekernel,
                              weights=newweights, varcov=varcov)
    raw <- divide.by.pixelarea(raw)
    smo <- raw
  }

  result <- if(is.im(smo)) smo[x$window, drop=FALSE]
            else solapply(smo, "[", i=x$window, drop=FALSE)

  # internal use only
  spill <- resolve.1.default(list(spill=FALSE), list(...))
  if(spill)
    return(list(result=result, sigma=sigma, varcov=varcov, raw = raw, edg=edg))

  # constrain values to be positive
  if(positive) 
    result <- posify(result)

  # normal return
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  attr(result, "kernel") <- kernel
  if(se)
    result <- list(estimate=result, SE=SE)
  return(result)
}

divideimage <- function(numer, denom) eval.im(numer/denom)

posify <- function(x, eps=.Machine$double.xmin) {
  force(eps) # scalpel
  if(is.im(x)) return(eval.im(pmax(eps, x)))
  if(inherits(x, "solist")) return(solapply(x, posify, eps=eps))
  if(is.numeric(x)) return(pmax(eps, x))
  # data frame or list
  if(is.list(x) && all(sapply(x, is.numeric)))
    return(lapply(x, posify, eps=eps))
  warning("Internal error: posify did not recognise data format")
  return(x)
}

divide.by.pixelarea <- function(x) {
  if(is.im(x)) {
    x$v <- x$v/(x$xstep * x$ystep)
  } else {
    for(i in seq_along(x))
      x[[i]]$v <- with(x[[i]], v/(xstep * ystep))
  }
  return(x)
}

denspppSEcalc <- function(x, sigma, varcov, ...,
                          weights, edge, diggle, at) {
  ## Calculate standard error, rather than estimate
  tau <- taumat <- NULL
  if(is.null(varcov)) {
    varconst <- 1/(4 * pi * prod(sigma))
    tau <- sigma/sqrt(2)
  } else {
    varconst <- 1/(4 * pi * sqrt(det(varcov)))
    taumat <- varcov/2
  }
  ## Calculate edge correction weights
  if(edge) {
    edgeim <- second.moment.calc(x, sigma, what="edge", ...,
                                 varcov=varcov)
    if(diggle || at == "points") {
      edgeX <- safelookup(edgeim, x, warn=FALSE)
      diggleX <- 1/edgeX
      diggleX[!is.finite(diggleX)] <- 0
    }
    edgeim <- edgeim[Window(x), drop=FALSE]
  }
  ## Perform smoothing
  if(!edge) {
    ## no edge correction
    V <- density(x, sigma=tau, varcov=taumat, ...,
                 weights=weights, edge=edge, diggle=diggle, at=at)
  } else if(!diggle) {
    ## edge correction e(u)
    V <- density(x, sigma=tau, varcov=taumat, ...,
                 weights=weights, edge=edge, diggle=diggle, at=at)
    V <- if(at == "pixels") (V/edgeim) else (V * diggleX)
  } else {
    ## Diggle edge correction e(x_i)
    wts <- diggleX * (weights %orifnull% 1)
    V <- density(x, sigma=tau, varcov=taumat, ...,
                 weights=wts, edge=edge, diggle=diggle, at=at)
  }
  V <- V * varconst
  return(sqrt(V))
}       


density.ppp

})

densitypointsEngine <- function(x, sigma, ...,
                                kernel="gaussian", 
                                scalekernel=is.character(kernel),
                                weights=NULL, edge=TRUE, varcov=NULL,
                                leaveoneout=TRUE, diggle=FALSE,
                                sorted=FALSE, spill=FALSE, cutoff=NULL) {
  debugging <- spatstat.options("developer")
  stopifnot(is.logical(leaveoneout))

  validate2Dkernel(kernel)
  if(is.character(kernel)) kernel <- match2DkernelName(kernel)
  isgauss <- identical(kernel, "gaussian")

  # constant factor in density computations
  if(is.null(varcov)) {
    const <- 1/sigma^2 
  } else {
    detSigma <- det(varcov)
    Sinv <- solve(varcov)
    const <- 1/sqrt(detSigma)
  }
  if(isgauss) {
    # absorb leading constant in Gaussian density
    const <- const/(2 * pi)
  }
  
  if(length(weights) == 0 || (!is.null(dim(weights)) && nrow(weights) == 0))
    weights <- NULL
  # Leave-one-out computation
  # cutoff: contributions from pairs of distinct points
  # closer than 8 standard deviations
  sd <- if(is.null(varcov)) sigma else sqrt(sum(diag(varcov)))
  if(is.null(cutoff)) 
    cutoff <- 8 * sd
  if(debugging)
    cat(paste("cutoff=", cutoff, "\n"))

  if(leaveoneout && npoints(x) > 1) {
    # ensure each point has its closest neighbours within the cutoff
    nndmax <- maxnndist(x)
    cutoff <- max(2 * nndmax, cutoff)
    if(debugging)
      cat(paste("adjusted cutoff=", cutoff, "\n"))
  }
  # validate weights
  if(is.null(weights)) {
    k <- 1
  } else if(is.matrix(weights) || is.data.frame(weights)) {
    k <- ncol(weights)
    stopifnot(nrow(weights) == npoints(x))
    weights <- as.data.frame(weights)
    weightnames <- colnames(weights)
  } else {
    k <- 1
    stopifnot(length(weights) == npoints(x) || length(weights) == 1)
  }
  # evaluate edge correction weights at points 
  if(edge) {
    win <- x$window
    if(isgauss && is.null(varcov) && win$type == "rectangle") {
      # evaluate Gaussian probabilities directly
      xr <- win$xrange
      yr <- win$yrange
      xx <- x$x
      yy <- x$y
      xprob <-
        pnorm(xr[2], mean=xx, sd=sigma) - pnorm(xr[1], mean=xx, sd=sigma)
      yprob <-
        pnorm(yr[2], mean=yy, sd=sigma) - pnorm(yr[1], mean=yy, sd=sigma)
      edgeweight <- xprob * yprob
    } else {
      edg <- second.moment.calc(x, sigma=sigma,
                                kernel=kernel,
                                scalekernel=scalekernel,
                                what="edge", varcov=varcov)
      edgeweight <- safelookup(edg, x, warn=FALSE)
    }
    if(diggle) {
      # Diggle edge correction
      # edgeweight is attached to each point
      if(is.null(weights)) {
        k <- 1
        weights <- 1/edgeweight
      } else {
        weights <- weights/edgeweight
      }
    }
  }

  if(isgauss &&
     spatstat.options("densityTransform") && spatstat.options("densityC")) {
    ## .................. experimental C code .....................
    if(debugging)
      cat('Using experimental code!\n')
    npts <- npoints(x)
    result <- if(k == 1) numeric(npts) else matrix(, npts, k)
    xx <- x$x
    yy <- x$y
    ## transform to standard coordinates
    if(is.null(varcov)) {
      xx <- xx/(sqrt(2) * sigma)
      yy <- yy/(sqrt(2) * sigma)
    } else {
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
    }
    if(is.null(weights)) {
      zz <- .C("Gdenspt",
               nxy     = as.integer(npts),
               x       = as.double(xx),
               y       = as.double(yy),
               rmaxi   = as.double(cutoff),
               result  = as.double(double(npts)))
      if(sorted) result <- zz$result else result[oo] <- zz$result
      result <- result * const
    } else if(k == 1) {
      wtsort <- if(sorted) weights else weights[oo]
      zz <- .C("Gwtdenspt",
               nxy     = as.integer(npts),
               x       = as.double(xx),
               y       = as.double(yy),
               rmaxi   = as.double(cutoff),
               weight  = as.double(wtsort),
               result  = as.double(double(npts)))
      if(sorted) result <- zz$result else result[oo] <- zz$result 
      result <- result * const
    } else {
      ## matrix of weights
      wtsort <- if(sorted) weights else weights[oo, ]
      for(j in 1:k) {
        zz <- .C("Gwtdenspt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 rmaxi   = as.double(cutoff),
                 weight  = as.double(wtsort[,j]),
                 result  = as.double(double(npts)))
        if(sorted) result[,j] <- zz$result else result[oo,j] <- zz$result
      }
      result <- result * const
    }
  } else if(isgauss && spatstat.options("densityC")) {
    # .................. C code ...........................
    if(debugging)
      cat('Using standard code.\n')
    npts <- npoints(x)
    result <- if(k == 1) numeric(npts) else matrix(, npts, k)
    # sort into increasing order of x coordinate (required by C code)
    if(sorted) {
      xx <- x$x
      yy <- x$y
    } else {
      oo <- fave.order(x$x)
      xx <- x$x[oo]
      yy <- x$y[oo]
    }
    if(is.null(varcov)) {
      # isotropic kernel
      if(is.null(weights)) {
        zz <- .C("denspt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 rmaxi   = as.double(cutoff),
                 sig     = as.double(sd),
                 result  = as.double(double(npts)))
        if(sorted) result <- zz$result else result[oo] <- zz$result 
      } else if(k == 1) {
        wtsort <- if(sorted) weights else weights[oo]
        zz <- .C("wtdenspt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 rmaxi   = as.double(cutoff),
                 sig     = as.double(sd),
                 weight  = as.double(wtsort),
                 result  = as.double(double(npts)))
        if(sorted) result <- zz$result else result[oo] <- zz$result 
       } else {
        # matrix of weights
        wtsort <- if(sorted) weights else weights[oo, ]
        for(j in 1:k) {
          zz <- .C("wtdenspt",
                   nxy     = as.integer(npts),
                   x       = as.double(xx),
                   y       = as.double(yy),
                   rmaxi   = as.double(cutoff),
                   sig     = as.double(sd),
                   weight  = as.double(wtsort[,j]),
                   result  = as.double(double(npts)))
          if(sorted) result[,j] <- zz$result else result[oo,j] <- zz$result
        }
      }
    } else {
      # anisotropic kernel
      flatSinv <- as.vector(t(Sinv))
      if(is.null(weights)) {
        zz <- .C("adenspt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 rmaxi   = as.double(cutoff),
                 detsigma = as.double(detSigma),
                 sinv    = as.double(flatSinv),
                 result  = as.double(double(npts)))
        if(sorted) result <- zz$result else result[oo] <- zz$result 
      } else if(k == 1) {
        # vector of weights
        wtsort <- if(sorted) weights else weights[oo]
        zz <- .C("awtdenspt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 rmaxi   = as.double(cutoff),
                 detsigma = as.double(detSigma),
                 sinv    = as.double(flatSinv),
                 weight  = as.double(wtsort),
                 result   = as.double(double(npts)))
        if(sorted) result <- zz$result else result[oo] <- zz$result 
      } else {
        # matrix of weights
        wtsort <- if(sorted) weights else weights[oo, ]
        for(j in 1:k) {
          zz <- .C("awtdenspt",
                   nxy     = as.integer(npts),
                   x       = as.double(xx),
                   y       = as.double(yy),
                   rmaxi   = as.double(cutoff),
                   detsigma = as.double(detSigma),
                   sinv    = as.double(flatSinv),
                   weight  = as.double(wtsort[,j]),
                   result  = as.double(double(npts)))
          if(sorted) result[,j] <- zz$result else result[oo,j] <- zz$result 
        }
      }
    }
  } else {
    # ..... interpreted code .........................................
    close <- closepairs(x, cutoff)
    i <- close$i
    j <- close$j
    d <- close$d
    npts <- npoints(x)
    result <- if(k == 1) numeric(npts) else matrix(, npts, k)
    # evaluate contribution from each close pair (i,j)
    if(isgauss) { 
      if(is.null(varcov)) {
        contrib <- const * exp(-d^2/(2 * sigma^2))
      } else {
        ## anisotropic kernel
        dx <- close$dx
        dy <- close$dy
        contrib <- const * exp(-(dx * (dx * Sinv[1,1] + dy * Sinv[1,2])
                                 + dy * (dx * Sinv[2,1] + dy * Sinv[2,2]))/2)
      }
    } else {
      contrib <- evaluate2Dkernel(kernel, close$dx, close$dy,
                                  sigma=sigma, varcov=varcov, ...)
    }
    ## sum (weighted) contributions
    ifac <- factor(i, levels=1:npts)
    if(is.null(weights)) {
      result <- tapply(contrib, ifac, sum)
    } else if(k == 1) {
      wcontrib <- contrib * weights[j]
      result <- tapply(wcontrib, ifac, sum)
    } else {
      for(kk in 1:k) {
        wcontribkk <- contrib * weights[j, kk]
        result[,kk] <- tapply(wcontribkk, ifac, sum)
      }
    }
    result[is.na(result)] <- 0
    #
  }
  # ----- contribution from point itself ----------------
  if(!leaveoneout) {
    # add contribution from point itself
    self <- const
    if(!is.null(weights))
      self <- self * weights
    result <- result + self
  }
  # ........  Edge correction ........................................
  if(edge && !diggle) 
    result <- result/edgeweight

  # ............. validate .................................
  npts <- npoints(x)
  if(k == 1) {
    result <- as.numeric(result)
    if(length(result) != npts) 
      stop(paste("Internal error: incorrect number of lambda values",
                 "in leave-one-out method:",
                 "length(lambda) = ", length(result),
                 "!=", npts, "= npoints"))
    if(anyNA(result)) {
      nwrong <- sum(is.na(result))
      stop(paste("Internal error:", nwrong, "NA or NaN",
                 ngettext(nwrong, "value", "values"),
                 "generated in leave-one-out method"))
    }
  } else {
    if(ncol(result) != k)
      stop(paste("Internal error: incorrect number of columns returned:",
                 ncol(result), "!=", k))
    colnames(result) <- weightnames
    if(nrow(result) != npts) 
      stop(paste("Internal error: incorrect number of rows of lambda values",
                 "in leave-one-out method:",
                 "nrow(lambda) = ", nrow(result),
                 "!=", npts, "= npoints"))
    if(anyNA(result)) {
      nwrong <- sum(!complete.cases(result))
      stop(paste("Internal error:", nwrong,
                 ngettext(nwrong, "row", "rows"),
                 "of NA values generated in leave-one-out method"))
    }
  }
  if(spill)
      return(list(result=result, sigma=sigma, varcov=varcov,
                  edg=edgeweight))
  # tack on bandwidth
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  # 
  return(result)
}

resolve.2D.kernel <- function(..., sigma=NULL, varcov=NULL, x, mindist=NULL,
                              adjust=1, bwfun=NULL, allow.zero=FALSE) {
  if(is.function(sigma)) {
    bwfun <- sigma
    sigma <- NULL
  }
  if(is.null(sigma) && is.null(varcov) && !is.null(bwfun)) {
    # call bandwidth selection function
    bw <- do.call.matched(bwfun, resolve.defaults(list(X=x), list(...)))
    # interpret the result as either sigma or varcov
    if(!is.numeric(bw))
      stop("bandwidth selector returned a non-numeric result")
    if(length(bw) %in% c(1,2)) {
      sigma <- as.numeric(bw)
      if(!all(sigma > 0)) {
        gripe <- "bandwidth selector returned negative value(s)"
        if(allow.zero) warning(gripe) else stop(gripe)
      }
    } else if(is.matrix(bw) && nrow(bw) == 2 && ncol(bw) == 2) {
      varcov <- bw
      if(!all(eigen(varcov)$values > 0))
        stop("bandwidth selector returned matrix with negative eigenvalues")
    } else stop("bandwidth selector did not return a matrix or numeric value")
  }
  sigma.given <- !is.null(sigma)
  varcov.given <- !is.null(varcov)
  if(sigma.given) {
    stopifnot(is.numeric(sigma))
    stopifnot(length(sigma) %in% c(1,2))
    if(!allow.zero)
      stopifnot(all(sigma > 0))
  }
  if(varcov.given)
    stopifnot(is.matrix(varcov) && nrow(varcov) == 2 && ncol(varcov)==2 )
  # reconcile
  ngiven <- varcov.given + sigma.given
  switch(ngiven+1,
         {
           # default
           w <- x$window
           sigma <- (1/8) * shortside(as.rectangle(w))
         },
         {
           if(sigma.given && length(sigma) == 2) 
             varcov <- diag(sigma^2)
           if(!is.null(varcov))
             sigma <- NULL
         },
         {
           stop(paste("Give only one of the arguments",
                      sQuote("sigma"), "and", sQuote("varcov")))
         })
  # apply adjustments
  if(!is.null(sigma))  sigma <- adjust * sigma
  if(!is.null(varcov)) varcov <- (adjust^2) * varcov
  #
  sd <- if(is.null(varcov)) sigma else sqrt(sum(diag(varcov)))
  cutoff <- 8 * sd
  uhoh <- if(!is.null(mindist) && cutoff < mindist) "underflow" else NULL
  result <- list(sigma=sigma, varcov=varcov, cutoff=cutoff, warnings=uhoh)
  return(result)
}


densitycrossEngine <- function(Xdata, Xquery, sigma, ...,
                               weights=NULL, edge=TRUE, varcov=NULL,
                               diggle=FALSE,
                               sorted=FALSE) {
  if(!is.null(varcov)) {
    detSigma <- det(varcov)
    Sinv <- solve(varcov)
  }
  if(length(weights) == 0 || (!is.null(dim(weights)) && nrow(weights) == 0))
    weights <- NULL
  ## Leave-one-out computation
  ## cutoff: contributions from pairs of distinct points
  ## closer than 8 standard deviations
  sd <- if(is.null(varcov)) sigma else sqrt(sum(diag(varcov)))
  cutoff <- 8 * sd
  # validate weights
  if(is.null(weights)) {
    k <- 1
  } else if(is.matrix(weights) || is.data.frame(weights)) {
    k <- ncol(weights)
    stopifnot(nrow(weights) == npoints(Xdata))
    weights <- as.data.frame(weights)
    weightnames <- colnames(weights)
  } else {
    k <- 1
    stopifnot(length(weights) == npoints(Xdata) || length(weights) == 1)
  }
  # evaluate edge correction weights at points 
  if(edge) {
    win <- Xdata$window
    if(diggle) {
      ## edge correction weights are attached to data points
      xedge <- Xdata
    } else {
      ## edge correction weights are applied at query points
      xedge <- Xquery
      if(!all(inside.owin(Xquery, , win)))
        stop(paste("Edge correction is not possible:",
                   "some query points lie outside the data window"),
             call.=FALSE)
    }
    if(is.null(varcov) && win$type == "rectangle") {
        ## evaluate Gaussian probabilities directly
      xr <- win$xrange
      yr <- win$yrange
      xx <- xedge$x
      yy <- xedge$y
      xprob <-
        pnorm(xr[2], mean=xx, sd=sigma) - pnorm(xr[1], mean=xx, sd=sigma)
      yprob <-
        pnorm(yr[2], mean=yy, sd=sigma) - pnorm(yr[1], mean=yy, sd=sigma)
      edgeweight <- xprob * yprob
    } else {
      edg <- second.moment.calc(Xdata, sigma=sigma,
                                what="edge", varcov=varcov)
      edgeweight <- safelookup(edg, xedge, warn=FALSE)
    }
    if(diggle) {
      ## Diggle edge correction
      ## edgeweight is attached to each data point
      if(is.null(weights)) {
        k <- 1
        weights <- 1/edgeweight
      } else {
        weights <- weights/edgeweight
      }
    }
  }
  
  ndata <- npoints(Xdata)
  nquery <- npoints(Xquery)
  result <- if(k == 1) numeric(nquery) else matrix(, nquery, k)
  ## coordinates
  xq <- Xquery$x
  yq <- Xquery$y
  xd <- Xdata$x
  yd <- Xdata$y
  if(!sorted) {
    ## sort into increasing order of x coordinate (required by C code)
    ooq <- fave.order(Xquery$x)
    xq <- xq[ooq]
    yq <- yq[ooq]
    ood <- fave.order(Xdata$x)
    xd <- xd[ood]
    yd <- yd[ood]
  }
  if(is.null(varcov)) {
    ## isotropic kernel
    if(is.null(weights)) {
      zz <- .C("crdenspt",
               nquery  = as.integer(nquery),
               xq      = as.double(xq),
               yq      = as.double(yq),
               ndata   = as.integer(ndata),
               xd      = as.double(xd),
               yd      = as.double(yd),
               rmaxi   = as.double(cutoff),
               sig     = as.double(sd),
               result  = as.double(double(nquery)))
      if(sorted) result <- zz$result else result[ooq] <- zz$result 
    } else if(k == 1) {
      wtsort <- if(sorted) weights else weights[ood]
      zz <- .C("wtcrdenspt",
               nquery  = as.integer(nquery),
               xq      = as.double(xq),
               yq      = as.double(yq),
               ndata   = as.integer(ndata),
               xd      = as.double(xd),
               yd      = as.double(yd),
               wd      = as.double(wtsort),
               rmaxi   = as.double(cutoff),
               sig     = as.double(sd),
               result  = as.double(double(nquery)))
      if(sorted) result <- zz$result else result[ooq] <- zz$result 
    } else {
      ## matrix of weights
      wtsort <- if(sorted) weights else weights[ood, ]
      for(j in 1:k) {
        zz <- .C("wtcrdenspt",
                 nquery  = as.integer(nquery),
                 xq      = as.double(xq),
                 yq      = as.double(yq),
                 ndata   = as.integer(ndata),
                 xd      = as.double(xd),
                 yd      = as.double(yd),
                 wd      = as.double(wtsort[,j]),
                 rmaxi   = as.double(cutoff),
                 sig     = as.double(sd),
                 result  = as.double(double(nquery)))
        if(sorted) result[,j] <- zz$result else result[ooq,j] <- zz$result
      }
      colnames(result) <- weightnames
    }
  } else {
    ## anisotropic kernel
    flatSinv <- as.vector(t(Sinv))
    if(is.null(weights)) {
      zz <- .C("acrdenspt",
               nquery  = as.integer(nquery),
               xq      = as.double(xq),
               yq      = as.double(yq),
               ndata   = as.integer(ndata),
               xd      = as.double(xd),
               yd      = as.double(yd),
               rmaxi   = as.double(cutoff),
               detsigma = as.double(detSigma),
               sinv    = as.double(flatSinv),
               result  = as.double(double(nquery)))
      if(sorted) result <- zz$result else result[ooq] <- zz$result 
    } else if(k == 1) {
      ## vector of weights
      wtsort <- if(sorted) weights else weights[ood]
      zz <- .C("awtcrdenspt",
               nquery  = as.integer(nquery),
               xq      = as.double(xq),
               yq      = as.double(yq),
               ndata   = as.integer(ndata),
               xd      = as.double(xd),
               yd      = as.double(yd),
               wd      = as.double(wtsort),
               rmaxi   = as.double(cutoff),
               detsigma = as.double(detSigma),
               sinv    = as.double(flatSinv),
               result   = as.double(double(nquery)))
      if(sorted) result <- zz$result else result[ooq] <- zz$result 
    } else {
      ## matrix of weights
      wtsort <- if(sorted) weights else weights[ood, ]
      for(j in 1:k) {
        zz <- .C("awtcrdenspt",
                 nquery  = as.integer(nquery),
                 xq      = as.double(xq),
                 yq      = as.double(yq),
                 ndata   = as.integer(ndata),
                 xd      = as.double(xd),
                 yd      = as.double(yd),
                 wd      = as.double(wtsort[,j]),
                 rmaxi   = as.double(cutoff),
                 detsigma = as.double(detSigma),
                 sinv    = as.double(flatSinv),
                 result  = as.double(double(nquery)))
        if(sorted) result[,j] <- zz$result else result[ooq,j] <- zz$result 
      }
      colnames(result) <- weightnames
    }
  }
  # ........  Edge correction ........................................
  if(edge && !diggle) 
    result <- result/edgeweight

  # tack on bandwidth
  attr(result, "sigma") <- sigma
  attr(result, "varcov") <- varcov
  # 
  return(result)
}


