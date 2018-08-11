#
# flipxy.R
#
# flip x and y coordinates
#
# $Revision: 1.3 $ $Date: 2017/02/07 07:22:47 $ 
#

flipxy <- function(X) {
  UseMethod("flipxy")
}

flipxy.ppp <- function(X) {
  stopifnot(is.ppp(X))
  ppp(X$y, X$x, marks=X$marks,
      window=flipxy(X$window), unitname=unitname(X),
      check=FALSE)
}

flipxypolygon <- function(p) {
  # flip x and y coordinates, and reinstate anticlockwise order
  oldy <- p$y
  p$y <- rev(p$x)
  p$x <- rev(oldy)
  # area and hole status unchanged
  return(p)
}

flipxy.owin <- function(X) {
  verifyclass(X, "owin")
  switch(X$type,
         rectangle={
           W <- owin(X$yrange, X$xrange, unitname=unitname(X))
         },
         polygonal={
           bdry <- lapply(X$bdry, flipxypolygon)
           W <- owin(poly=bdry, check=FALSE, unitname=unitname(X))
         },
         mask={
           W <- owin(mask=t(X$m),
                     xy=list(x=X$yrow, y=X$xcol),
                     unitname=unitname(X))
         },
         stop("Unrecognised window type")
         )
  return(W)
}

flipxy.psp <- function(X) {
  stopifnot(is.psp(X))
  flipends <- (X$ends)[, c(2L,1L,4L,3L), drop=FALSE]
  as.psp(flipends, window=flipxy(X$window), marks=X$marks,
         unitname=unitname(X), check=FALSE)
}

flipxy.im <- function(X) {
  im(t(X$v), xcol=X$yrow, yrow=X$xcol, unitname=unitname(X))
}

