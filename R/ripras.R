#
#	ripras.S	Ripley-Rasson estimator of domain
#
#
#	$Revision: 1.14 $	$Date: 2014/10/24 00:22:30 $
#
#
#
#
#-------------------------------------
bounding.box.xy <- function(x, y=NULL) {
  xy <- xy.coords(x,y)
  if(length(xy$x) == 0)
    return(NULL)
  owin(range(xy$x), range(xy$y), check=FALSE)
}

convexhull.xy <- function(x, y=NULL) {
  xy <- xy.coords(x, y)
  x <- xy$x
  y <- xy$y
  if(length(x) < 3)
    return(NULL)
  h <- rev(chull(x, y))  # must be anticlockwise
  if(length(h) < 3)
    return(NULL)
  w <- owin(poly=list(x=x[h], y=y[h]), check=FALSE)
  return(w)
}

ripras <- function(x, y=NULL, shape="convex", f) {
  xy <- xy.coords(x, y)
  n <- length(xy$x)
  w <- switch(shape,
              convex = convexhull.xy(xy),
              rectangle = boundingbox(xy),
              stop(paste("Unrecognised option: shape=", dQuote(shape))))
  if(is.null(w))
    return(NULL)
  # expansion factor
  if(!missing(f))
    stopifnot(is.numeric(f) && length(f) == 1 && f >= 1)
  else switch(shape,
              convex = {
                # number of vertices
                m <- summary(w)$nvertices
                f <- if(m < n) 1/sqrt(1 - m/n) else 2
              },
              rectangle = {
                f <- (n+1)/(n-1)
              })
  # centroid
  ce <- unlist(centroid.owin(w))
  # shift centroid to origin
  W <- shift(w, -ce)
  # rescale
  W <- affine(W, mat=diag(c(f,f)))
  # shift origin to centroid
  W <- shift(W, ce)
  return(W)
}

