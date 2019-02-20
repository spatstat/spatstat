#
#    discarea.R
#
#  $Revision: 1.21 $  $Date: 2019/02/20 03:34:50 $
#
#
#  Compute area of intersection between a disc and a window,
#
discpartarea <- function(X, r, W=as.owin(X)) {
  if(!missing(W)) {
    verifyclass(W, "owin")
    if(!inherits(X, "ppp"))
      X <- as.ppp(X, W)
  }
  verifyclass(X, "ppp")

  n <- X$n
  if(is.matrix(r) && nrow(r) != n)
    stop("the number of rows of r should match the number of points in X")
  if(!is.matrix(r)) {
    nr <- length(r)
    r <- matrix(r, nrow=n, ncol=nr, byrow=TRUE)
  } else {
    nr <- ncol(r)
  }
  
  W <- as.polygonal(W)
  
  # convert polygon to line segments
  Y <- edges(W)
  # remove vertical segments (contribution is zero)
  vert <- (Y$ends$x1 == Y$ends$x0)
  Y <- Y[!vert]
  # go
  z <- .C("discareapoly",
          nc=as.integer(n),
          xc=as.double(X$x),
          yc=as.double(X$y),
          nr=as.integer(nr),
          rmat=as.double(r),
          nseg=as.integer(Y$n),
          x0=as.double(Y$ends$x0),
          y0=as.double(Y$ends$y0),
          x1=as.double(Y$ends$x1),
          y1=as.double(Y$ends$y1),
          eps=as.double(.Machine$double.eps),
          out=as.double(numeric(length(r))),
          PACKAGE = "spatstat")
  areas <- matrix(z$out, n, nr)
  return(areas)
}

# Compute area of dilation of point pattern
# using Dirichlet tessellation or distmap
#  (areas of other dilations using distmap)

dilated.areas <- function(X, r, W=as.owin(X), ...,
                          constrained=TRUE,
                          exact=FALSE) {
  if(is.matrix(r)) {
    if(sum(dim(r) > 1) > 1L)
      stop("r should be a vector or single value")
    r <- as.vector(r)
  }
  if(exact && !is.ppp(X)) {
    exact <- FALSE
    warning("Option exact=TRUE is only available for ppp objects")
  }
  if(!constrained || is.null(W)) {
    # unconstrained dilation
    bb <- as.rectangle(X)
    W <- grow.rectangle(bb, max(r))
    if(is.owin(X))
      X <- rebound.owin(X, W)
    else
      X$window <- W
  } else W <- as.owin(W)
  if(!exact) {
    D <- distmap(X, ...)
    pixelarea <- D$xstep * D$ystep
    Dvals <- D[W, drop=TRUE]
    if(is.im(Dvals))
      Dvals <- as.vector(as.matrix(Dvals))
    Dvals <- Dvals[!is.na(Dvals)]
    rr <- c(-1, r)
    h <- cumsum(whist(Dvals, rr))
    return(h * pixelarea)
  }
  npts <- npoints(X)
  nr <- length(r)
  if(npts == 0)
    return(numeric(nr))
  else if(npts == 1L) 
    return(discpartarea(X, r, W))
  samebox <- (W$type == "rectangle") &&
              isTRUE(all.equal(W, as.owin(X)))
  needclip <- constrained && !samebox
  X <- unique(X)
  dd <- dirichlet(X)
  til <- tiles(dd)
  #' some data points may not have a tile
  whichpoint <- as.integer(names(til))
  partareas <- matrix(0, length(til), nr)
  for(j in seq_along(til)) {
    Tj <- til[[j]]
    if(needclip)
      Tj <- intersect.owin(Tj, W)
    i <- whichpoint[j]
    partareas[j,] <- discpartarea(X[i], r, Tj)
  }
  return(colSums(partareas))
}

  
