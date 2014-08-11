#
#
#      quadscheme.S
#
#      $Revision: 4.25 $    $Date: 2010/03/08 08:23:04 $
#
#      quadscheme()    generate a quadrature scheme from 
#		       data and dummy point patterns.
#
#      quadscheme.spatial()    case where both patterns are unmarked
#
#      quadscheme.replicated() case where data are multitype
#
#
#---------------------------------------------------------------------

quadscheme <- function(data, dummy, method="grid", ...) {
        #
	# generate a quadrature scheme from data and dummy patterns.
	#
	# Other arguments control how the quadrature weights are computed
        #

  data <- as.ppp(data)

  if(missing(dummy)) {
    # create dummy points
    dummy <- default.dummy(data, method=method, ...)
    # extract full set of parameters used to create dummy points
    dp <- attr(dummy, "dummy.parameters")
    # extract recommended parameters for computing weights
    wp <- attr(dummy, "weight.parameters")
  } else {
    # user-supplied dummy points
    if(!is.ppp(dummy)) {
      # convert to ppp object
      dummy <- as.ppp(dummy, data$window, check=FALSE)
      # confine dummy points to data window 
      dummy <- dummy[data$window]
      wp <- dp <- list()
    } else {
      # if it's already a ppp, it may have been created by default.dummy
     dp <- attr(dummy, "dummy.parameters")
     wp <- attr(dummy, "weight.parameters")     
   }
  }
  # arguments supplied directly to quadscheme()
  # override any arguments passed as attributes
  wp <- resolve.defaults(list(method=method), list(...), wp)
  
  mX <- is.marked(data)
  mD <- is.marked(dummy)

  if(!mX && !mD)
    Q <- do.call("quadscheme.spatial",
                 append(list(data, dummy, check=FALSE), wp))
  else if(mX && !mD)
    Q <- do.call("quadscheme.replicated",
                 append(list(data, dummy, check=FALSE), wp))
  else if(!mX && mD)
    stop("dummy points are marked but data are unmarked")
  else
    stop("marked data and marked dummy points -- sorry, this case is not implemented")

  # record parameters used to make dummy points
  Q$param$dummy <- dp

  return(Q)
}

quadscheme.spatial <-
  function(data, dummy, method="grid", ...) {
        #
	# generate a quadrature scheme from data and dummy patterns.
	#
	# The 'method' may be "grid" or "dirichlet"
	#
	# '...' are passed to gridweights() or dirichlet.weights()
        #
        # quadscheme.spatial:
        #       for unmarked point patterns.
        #
        #       weights are determined only by spatial locations
        #       (i.e. weight computations ignore any marks)
	#
        # No two points should have the same spatial location
        # 

    check <- resolve.defaults(list(...), list(check=TRUE))$check
    
	data <- as.ppp(data, check=check)
        dummy <- as.ppp(dummy, data$window, check=check)
		# note data$window is the DEFAULT quadrature window
		# applicable when 'dummy' does not contain a window

        if(is.marked(data, dfok=TRUE))
          warning("marks in data pattern - ignored")
        if(is.marked(dummy, dfok=TRUE))
          warning("marks in dummy pattern - ignored")
        
	both <- as.ppp(concatxy(data, dummy), dummy$window, check=check)
	switch(method,
		grid={
			w <- gridweights(both, window= dummy$window, ...)
		},
		dirichlet = {
			w <- dirichlet.weights(both, window=dummy$window, ...)
		},
		{ 
			stop(paste("unrecognised method", sQuote(method)))
		}
	)

        # parameters actually used to make weights
        wp <- attr(w, "weight.parameters")
        param <- list(weight = wp, dummy = NULL)

	Q <- quad(data, dummy, w, param)
        return(Q)
}

"quadscheme.replicated" <-
  function(data, dummy, method="grid", ...) {
        #
	# generate a quadrature scheme from data and dummy patterns.
	#
	# The 'method' may be "grid" or "dirichlet"
	#
	# '...' are passed to gridweights() or dirichlet.weights()
        #
        # quadscheme.replicated:
        #       for multitype point patterns.
        #
        # No two points in 'data'+'dummy' should have the same spatial location

    check <- resolve.defaults(list(...), list(check=TRUE))$check

    data <- as.ppp(data, check=check)
    dummy <- as.ppp(dummy, data$window, check=check)
		# note data$window is the DEFAULT quadrature window
		# unless otherwise specified in 'dummy'

        if(!is.marked(data))
          stop("data pattern does not have marks")
        if(is.marked(dummy, dfok=TRUE))
          warning("dummy points have marks --- ignored")

        # first, ignore marks and compute spatial weights
        P <- quadscheme.spatial(unmark(data), dummy, method, ...)
        W <- w.quad(P)
        iz <- is.data(P)
        Wdat <- W[iz]
        Wdum <- W[!iz]

        # find the set of all possible marks

        if(!is.multitype(data))
          stop("data pattern is not multitype")
        data.marks <- marks(data)
        markset <- levels(data.marks)
        nmarks <- length(markset)

        # replicate dummy points, one copy for each possible mark
        # -> dummy x {1,..,K}
        
        dumdum <- cartesian(dummy, markset)
        Wdumdum <- rep(Wdum, nmarks)
        
        # also make dummy marked points at same locations as data points
        # but with different marks

        dumdat <- cartesian(unmark(data), markset)
        Wdumdat <- rep(Wdat, nmarks)
        Mdumdat <- marks(dumdat)
        
        Mrepdat <- rep(data.marks, nmarks)

        ok <- (Mdumdat != Mrepdat)
        dumdat <- dumdat[ok,]
        Wdumdat <- Wdumdat[ok]

        # combine the two dummy patterns

        dumb <- superimpose(dumdum, dumdat, W=dummy$window)
        Wdumb <- c(Wdumdum, Wdumdat)

        # record the quadrature parameters
        param <- list(weight = P$param$weight, dummy = NULL)

        # wrap up
	Q <- quad(data, dumb, c(Wdat, Wdumb), param)
        return(Q)
}


"cartesian" <-
function(pp, markset, fac=TRUE) {
  # given an unmarked point pattern 'pp'
  # and a finite set of marks,
  # create the marked point pattern which is
  # the Cartesian product, consisting of all pairs (u,k)
  # where u is a point of 'pp' and k is a mark in 'markset'
  nmarks <- length(markset)
  result <- ppp(rep(pp$x, nmarks),
                rep(pp$y, nmarks),
                window=pp$window,
                check=FALSE)
  marx <- rep(markset, rep(pp$n, nmarks))
  if(fac)
    marx <- factor(marx, levels=markset)
  marks(result) <- marx
  return(result)
}


validate.quad <- function(Q, fatal=FALSE, repair=TRUE, announce=FALSE) {
  X <- Q$data
  D <- Q$dummy
  mX <- is.marked(X)
  mD <- is.marked(D)
  nbg <- function(whinge, fatal=FALSE, announce=FALSE) {
    if(fatal)
      stop(whinge, call.=FALSE)
    else {
      if(announce)
        warning(whinge, call.=FALSE)
      return(FALSE)
    }
  }
  if(mX != mD) {
    whinge <-
      if(mX)
        "data points are marked, but dummy points are not"
      else
        "dummy points are marked, but data points are not"
    return(nbg(whinge, fatal, announce))
  }
  if(!mX)
    return(TRUE)
  # marked points 
  fX <- is.factor(Xmarx <- marks(X))
  fD <- is.factor(Dmarx <- marks(D))
  if(fX != fD) {
    whinge <-
      if(fX)
        "data points are multitype, but dummy points are not"
      else
        "dummy points are multitype, but data points are not"
    return(nbg(whinge, fatal, announce))
  }
  if(!fX)
    return(TRUE)
  # multitype points
  lX <- levels(Xmarx)
  lD <- levels(Dmarx)
  if(length(lX) != length(lD) || any(lX != lD)) {
    whinge <- "data and dummy points have different sets of possible marks"
    return(nbg(whinge, fatal, announce))
  }
  return(TRUE)
}

  

pixelquad <- function(X, W=as.owin(X)) {
  # make a quadscheme with a dummy point at every pixel
  verifyclass(X, "ppp")

  # convert window to mask if not already one
  W <- as.owin(W)
  M <- as.mask(W)
  MM <- M$m
  pixelarea <- M$xstep * M$ystep
  
  # create pixel coordinates and corresponding row, column indices
  xx <- as.vector(raster.x(M)[MM])
  yy <- as.vector(raster.y(M)[MM])
  cc <- as.vector(col(MM)[MM])
  rr <- as.vector(row(MM)[MM])
  Nr <- M$dim[1]
  Nc <- M$dim[2]
  
  # discretise data points
  ij <- nearest.raster.point(X$x, X$y, M)
  ijrow <- ij$row
  ijcol <- ij$col

  # tabulate pixel locations of data points
  Xtab <- table(row=factor(ijrow, levels=1:Nr),
                col=factor(ijcol, levels=1:Nc))

  # every pixel contains exactly one dummy point,
  # so the total count of quadrature points in each pixel is:
  Qtab <- Xtab + 1

  # compute counting weights for data points
  wdat <- 1/Qtab[cbind(ijrow, ijcol)]
  # compute counting weights for dummy points
  wdum <- 1/Qtab[cbind(rr, cc)]

  wboth <- pixelarea * c(wdat, wdum)
  
  # create quadrature scheme
  dum <- ppp(xx, yy, window=W, check=FALSE)
  Q <- quad(X, dum, wboth)

  attr(Q, "M") <- M
  return(Q)
}

  
  
