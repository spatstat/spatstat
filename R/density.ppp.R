#
#  density.ppp.R
#
#  Method for 'density' for point patterns
#
#  $Revision: 1.110 $    $Date: 2020/05/09 04:02:00 $
#

# ksmooth.ppp <- function(x, sigma, ..., edge=TRUE) {
#  .Deprecated("density.ppp", package="spatstat")
#  density.ppp(x, sigma, ..., edge=edge)
# }

density.ppp <- local({
  
density.ppp <- function(x, sigma=NULL, ...,
                        weights=NULL, edge=TRUE, varcov=NULL,
                        at="pixels", leaveoneout=TRUE,
                        adjust=1, diggle=FALSE, se=FALSE, 
                        kernel="gaussian",
                        scalekernel=is.character(kernel),
                        positive=FALSE, verbose=TRUE) {
  verifyclass(x, "ppp")

  output <- pickoption("output location type", at,
                       c(pixels="pixels",
                         points="points"))

  if(!identical(kernel, "gaussian")) {
    validate2Dkernel(kernel)
    ## kernel is only partly implemented!
    if(se)
      stop("Standard errors are not implemented for non-Gaussian kernel")
    if(verbose && scalekernel &&
       (is.function(sigma) || (is.null(sigma) && is.null(varcov))))
      warning("Bandwidth selection will be based on Gaussian kernel")
  }
  
  ker <- resolve.2D.kernel(..., sigma=sigma, varcov=varcov, x=x, adjust=adjust)
  sigma <- ker$sigma
  varcov <- ker$varcov
  ## sigma.is.infinite <- ker$infinite
  
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

  ## infinite bandwidth
  if(bandwidth.is.infinite(sigma)) {
    #' uniform estimate
    nx <- npoints(x)
    single <- is.null(dim(weights))
    totwt <- if(is.null(weights)) nx else
             if(single) sum(weights) else colSums(weights)
    if(!edge) totwt <- 0 * totwt
    W <- Window(x)
    A <- area.owin(W)
    switch(output,
           pixels = {
             E <- solapply(totwt/A, as.im, W=W, ...)
             names(E) <- colnames(weights)
             if(single) E <- E[[1L]]
           },
           points = {
             numerator <- rep(totwt, each=nx)
             if(!single) numerator <- matrix(numerator, nrow=nx)
             if(leaveoneout && edge) 
               numerator <- numerator - (weights %orifnull% 1)
             E <- numerator/A
             if(!single)
               colnames(E) <- colnames(weights)
           })
    result <- if(se) list(estimate=E, SE=SE) else E
    return(result)
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
    if(verbose && !is.null(uhoh <- attr(result, "warnings"))) {
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
    ## Math.im / Math.imlist not yet working
    smo <- imagelistOp(raw, edg, "/")
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
  nx <- npoints(x)

  single <- is.null(dim(weights))
  
  if(bandwidth.is.infinite(sigma)) {
    #' special case - uniform
    totwt2 <- if(is.null(weights)) nx else
              if(single) sum(weights^2) else colSums(weights^2)
    if(!edge)
      totwt2 <- 0 * totwt2
    W <- Window(x)
    A <- area.owin(W)
    switch(at,
           pixels = {
             V <- solapply(totwt2/A, as.im, W=W, ...)
             names(V) <- colnames(weights)
             if(single) V <- V[[1L]]
           },
           points = {
             numerator <- rep(totwt2, each=nx)
             if(!single) numerator <- matrix(numerator, nrow=nx)
             leaveoneout <- resolve.1.default(list(leaveoneout=TRUE), list(...))
             if(edge && leaveoneout) 
               numerator <- numerator - (weights %orifnull% 1)^2
             V <- numerator/A
             if(!single) 
               colnames(V) <- colnames(weights)
           })
    return(sqrt(V))
  }

  ## Usual case
  tau <- taumat <- NULL
  if(is.null(varcov)) {
    varconst <- 1/(4 * pi * prod(ensure2vector(sigma)))
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
    V <- if(at == "points") V * diggleX else imagelistOp(V, edgeim, "/")
  } else {
    ## Diggle edge correction e(x_i)
    wts <- if(is.null(weights)) diggleX else (diggleX * weights)
    V <- density(x, sigma=tau, varcov=taumat, ...,
                 weights=wts, edge=edge, diggle=diggle, at=at)
  }
  V <- V * varconst
  return(sqrt(V))
}       


density.ppp

})

densitypointsEngine <- function(x, sigma=NULL, ...,
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

  if(isgauss) {
    ## constant factor in Gaussian density
    if(is.null(varcov)) {
      gaussconst <- 1/(2 * pi * sigma^2)
    } else {
      detSigma <- det(varcov)
      Sinv <- solve(varcov)
      gaussconst <- 1/(2 * pi * sqrt(detSigma))
    }
  }
  
  if(length(weights) == 0 || (!is.null(dim(weights)) && nrow(weights) == 0))
    weights <- NULL
  
  ## infinite bandwidth
  if(bandwidth.is.infinite(sigma)) {
    #' uniform estimate
    nx <- npoints(x)
    single <- is.null(dim(weights))
    totwt <- if(is.null(weights)) nx else
             if(single) sum(weights) else colSums(weights)
    if(!edge) totwt <- 0 * totwt
    W <- Window(x)
    A <- area.owin(W)
    numerator <- rep(totwt, each=nx)
    if(!single) numerator <- matrix(numerator, nrow=nx)
    if(leaveoneout && edge) 
      numerator <- numerator - (weights %orifnull% 1)
    result <- numerator/A
    if(!single)
      colnames(result) <- colnames(weights)
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

  if(leaveoneout && npoints(x) > 1) {
    ## ensure each point has its closest neighbours within the cutoff
    nndmax <- maxnndist(x)
    cutoff <- max(2 * nndmax, cutoff)
    if(debugging)
      cat(paste("adjusted cutoff=", cutoff, "\n"))
  }
  # validate weights
  if(is.null(weights)) {
    k <- 1L
  } else if(is.matrix(weights) || is.data.frame(weights)) {
    k <- ncol(weights)
    stopifnot(nrow(weights) == npoints(x))
    weights <- as.data.frame(weights)
    weightnames <- colnames(weights)
  } else {
    k <- 1L
    stopifnot(length(weights) == npoints(x) || length(weights) == 1L)
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
        pnorm(xr[2L], mean=xx, sd=sigma) - pnorm(xr[1L], mean=xx, sd=sigma)
      yprob <-
        pnorm(yr[2L], mean=yy, sd=sigma) - pnorm(yr[1L], mean=yy, sd=sigma)
      edgeweight <- xprob * yprob
    } else {
      edg <- second.moment.calc(x, sigma=sigma,
                                kernel=kernel,
                                scalekernel=scalekernel,
                                what="edge", varcov=varcov, ...)
      edgeweight <- safelookup(edg, x, warn=FALSE)
    }
    if(diggle) {
      # Diggle edge correction
      # edgeweight is attached to each point
      if(is.null(weights)) {
        k <- 1L
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
    result <- if(k == 1L) numeric(npts) else matrix(, npts, k)
    xx <- x$x
    yy <- x$y
    ## transform to standard coordinates
    if(is.null(varcov)) {
      xx <- xx/(sqrt(2) * sigma)
      yy <- yy/(sqrt(2) * sigma)
    } else {
      xy <- cbind(xx, yy) %*% matrixsqrt(Sinv/2)
      xx <- xy[,1L]
      yy <- xy[,2L]
      sorted <- FALSE
    }
    ## cutoff in standard coordinates
    sd <- sigma %orifnull% sqrt(min(eigen(varcov)$values))
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
               result  = as.double(double(npts)),
               PACKAGE = "spatstat")
      if(sorted) result <- zz$result else result[oo] <- zz$result
      result <- result * gaussconst
    } else if(k == 1L) {
      wtsort <- if(sorted) weights else weights[oo]
      zz <- .C("Gwtdenspt",
               nxy     = as.integer(npts),
               x       = as.double(xx),
               y       = as.double(yy),
               rmaxi   = as.double(cutoff),
               weight  = as.double(wtsort),
               result  = as.double(double(npts)),
               PACKAGE = "spatstat")
      if(sorted) result <- zz$result else result[oo] <- zz$result 
      result <- result * gaussconst
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
                 result  = as.double(double(npts)),
                 PACKAGE = "spatstat")
        if(sorted) result[,j] <- zz$result else result[oo,j] <- zz$result
      }
      result <- result * gaussconst
    }
  } else if(isgauss && spatstat.options("densityC")) {
    # .................. C code ...........................
    if(debugging)
      cat('Using standard code.\n')
    npts <- npoints(x)
    result <- if(k == 1L) numeric(npts) else matrix(, npts, k)
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
                 sig     = as.double(sigma),
                 result  = as.double(double(npts)),
                 PACKAGE = "spatstat")
        if(sorted) result <- zz$result else result[oo] <- zz$result 
      } else if(k == 1L) {
        wtsort <- if(sorted) weights else weights[oo]
        zz <- .C("wtdenspt",
                 nxy     = as.integer(npts),
                 x       = as.double(xx),
                 y       = as.double(yy),
                 rmaxi   = as.double(cutoff),
                 sig     = as.double(sigma),
                 weight  = as.double(wtsort),
                 result  = as.double(double(npts)),
                 PACKAGE = "spatstat")
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
                   sig     = as.double(sigma),
                   weight  = as.double(wtsort[,j]),
                   result  = as.double(double(npts)),
                   PACKAGE = "spatstat")
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
                 result  = as.double(double(npts)),
                 PACKAGE = "spatstat")
        if(sorted) result <- zz$result else result[oo] <- zz$result 
      } else if(k == 1L) {
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
                 result   = as.double(double(npts)),
                 PACKAGE = "spatstat")
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
                   result  = as.double(double(npts)),
                   PACKAGE = "spatstat")
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
    result <- if(k == 1L) numeric(npts) else matrix(, npts, k)
    # evaluate contribution from each close pair (i,j)
    if(isgauss) { 
      if(is.null(varcov)) {
        contrib <- gaussconst * exp(-d^2/(2 * sigma^2))
      } else {
        ## anisotropic kernel
        dx <- close$dx
        dy <- close$dy
        contrib <- gaussconst * exp(-(dx * (dx * Sinv[1L,1L] + dy * Sinv[1L,2L])
                               + dy * (dx * Sinv[2L,1L] + dy * Sinv[2L,2L]))/2)
      }
    } else {
      contrib <- evaluate2Dkernel(kernel, close$dx, close$dy,
                                  sigma=sigma, varcov=varcov,
                                  scalekernel=scalekernel, ...)
    }
    ## sum (weighted) contributions
    ## query point i, data point j
    ifac <- factor(i, levels=1:npts)
    if(is.null(weights)) {
      result <- tapplysum(contrib, list(ifac))
    } else if(k == 1L) {
      wcontrib <- contrib * weights[j]
      result <- tapplysum(wcontrib, list(ifac))
    } else {
      for(kk in 1:k) {
        wcontribkk <- contrib * weights[j, kk]
        result[,kk] <- tapplysum(wcontribkk, list(ifac))
      }
    }
    #
  }
  # ----- contribution from point itself ----------------
  if(!leaveoneout) {
    #' add contribution from point itself
    if(isgauss) {
      self <- gaussconst
    } else {
      self <- evaluate2Dkernel(kernel, 0, 0, sigma=sigma, varcov=varcov,
                               scalekernel=scalekernel, ...)
    }
    if(!is.null(weights))
      self <- self * weights
    result <- result + self
  }
  # ........  Edge correction ........................................
  if(edge && !diggle) 
    result <- result/edgeweight

  # ............. validate .................................
  npts <- npoints(x)
  if(k == 1L) {
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
    if(length(bw) %in% c(1L,2L)) {
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
    stopifnot(length(sigma) %in% c(1L,2L))
    if(!allow.zero)
      stopifnot(all(sigma > 0))
  }
  if(varcov.given)
    stopifnot(is.matrix(varcov) && nrow(varcov) == 2 && ncol(varcov)==2 )
  # reconcile
  ngiven <- varcov.given + sigma.given
  switch(ngiven+1L,
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

densitycrossEngine <- function(Xdata, Xquery, sigma=NULL, ...,
                               kernel="gaussian", 
                               scalekernel=is.character(kernel),
                               weights=NULL, edge=TRUE, varcov=NULL,
                               diggle=FALSE,
                               sorted=FALSE, cutoff=NULL) {
  validate2Dkernel(kernel)
  if(is.character(kernel)) kernel <- match2DkernelName(kernel)
  isgauss <- identical(kernel, "gaussian") && scalekernel

  if(length(weights) == 0 || (!is.null(dim(weights)) && nrow(weights) == 0))
    weights <- NULL
  # validate weights
  if(is.null(weights)) {
    k <- 1L
  } else if(is.matrix(weights) || is.data.frame(weights)) {
    k <- ncol(weights)
    stopifnot(nrow(weights) == npoints(Xdata))
    weights <- as.data.frame(weights)
    weightnames <- colnames(weights)
  } else {
    k <- 1L
    stopifnot(length(weights) == npoints(Xdata) || length(weights) == 1L)
  }

  #' infinite bandwidth
  if(bandwidth.is.infinite(sigma)) {
    #' uniform estimate
    single <- is.null(dim(weights))
    totwt <- if(is.null(weights)) npoints(Xdata) else
             if(single) sum(weights) else colSums(weights)
    if(!edge) totwt <- 0 * totwt
    lam <- totwt/area.owin(Window(Xdata))
    result <- if(single) rep(lam, npoints(Xquery)) else 
              matrix(lam, npoints(Xquery), length(lam), byrow=TRUE,
                     dimnames=list(NULL, colnames(weights)))
    return(result)
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
    if(isgauss && is.null(varcov) && win$type == "rectangle") {
      ## evaluate Gaussian probabilities directly
      xr <- win$xrange
      yr <- win$yrange
      xx <- xedge$x
      yy <- xedge$y
      xprob <-
        pnorm(xr[2L], mean=xx, sd=sigma) - pnorm(xr[1L], mean=xx, sd=sigma)
      yprob <-
        pnorm(yr[2L], mean=yy, sd=sigma) - pnorm(yr[1L], mean=yy, sd=sigma)
      edgeweight <- xprob * yprob
    } else {
      ## evaluate edge correction weights
      edg <- second.moment.calc(Xdata, what="edge",
                                kernel=kernel, scalekernel=scalekernel,
                                sigma=sigma, varcov=varcov)
      edgeweight <- safelookup(edg, xedge, warn=FALSE)
    }
    if(diggle) {
      ## Diggle edge correction
      ## edgeweight is attached to each data point
      if(is.null(weights)) {
        k <- 1L
        weights <- 1/edgeweight
      } else {
        weights <- weights/edgeweight
      }
    }
  }

  ## cutoff distance (beyond which the kernel value is treated as zero)
  ## NB: input argument 'cutoff' is either NULL or
  ##     an absolute distance (if scalekernel=FALSE)
  ##     a number of standard deviations (if scalekernel=TRUE)
  cutoff <- cutoff2Dkernel(kernel, sigma=sigma, varcov=varcov,
                           scalekernel=scalekernel, cutoff=cutoff,
                           fatal=TRUE)
  ## cutoff is now an absolute distance
  
  ndata <- npoints(Xdata)
  nquery <- npoints(Xquery)
  if(!isgauss) {
    ## .................. non-Gaussian kernel ........................
    close <- crosspairs(Xdata, Xquery, cutoff)
    contrib <- evaluate2Dkernel(kernel, close$dx, close$dy,
                                sigma=sigma, varcov=varcov,
                                scalekernel=scalekernel, ...)
    ## sum the (weighted) contributions
    i <- close$i
    j <- close$j
    jfac <- factor(j, levels=seq_len(nquery))
    if(is.null(weights)) {
      result <- tapplysum(contrib, list(jfac))
    } else if(k == 1L) {
      wcontrib <- contrib * weights[i]
      result <- tapplysum(wcontrib, list(jfac))
    } else {
      result <- matrix(, nquery, k)
      for(kk in 1:k) {
        wcontribkk <- contrib * weights[i, kk]
        result[,kk] <- tapplysum(wcontribkk, list(jfac))
      }
    }
  } else {
    ## ................. Gaussian kernel ...................
    result <- if(k == 1L) numeric(nquery) else matrix(, nquery, k)
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
                 sig     = as.double(sigma),
                 result  = as.double(double(nquery)),
                 PACKAGE = "spatstat")
        if(sorted) result <- zz$result else result[ooq] <- zz$result 
      } else if(k == 1L) {
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
                 sig     = as.double(sigma),
                 result  = as.double(double(nquery)),
                 PACKAGE = "spatstat")
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
                   sig     = as.double(sigma),
                   result  = as.double(double(nquery)),
                   PACKAGE = "spatstat")
          if(sorted) result[,j] <- zz$result else result[ooq,j] <- zz$result
        }
        colnames(result) <- weightnames
      }
    } else {
      ## anisotropic kernel
      detSigma <- det(varcov)
      Sinv <- solve(varcov)
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
                 result  = as.double(double(nquery)),
                 PACKAGE = "spatstat")
        if(sorted) result <- zz$result else result[ooq] <- zz$result 
      } else if(k == 1L) {
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
                 result   = as.double(double(nquery)),
                 PACKAGE = "spatstat")
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
                   result  = as.double(double(nquery)),
                   PACKAGE = "spatstat")
          if(sorted) result[,j] <- zz$result else result[ooq,j] <- zz$result 
        }
        colnames(result) <- weightnames
      }
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

bandwidth.is.infinite <- function(sigma) {
  sigma <- as.numeric(sigma)
  return((length(sigma) > 0) && all(sigma == Inf))
}
