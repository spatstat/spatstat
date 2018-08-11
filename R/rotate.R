#
#	rotate.S
#
#	$Revision: 1.21 $	$Date: 2014/10/24 00:22:30 $
#

rotxy <- function(X, angle=pi/2) {
  co <- cos(angle)
  si <- sin(angle)
  list(x = co * X$x - si * X$y,
       y = si * X$x + co * X$y)
}

rotxypolygon <- function(p, angle=pi/2) {
  p[c("x","y")] <- rotxy(p, angle=angle)
  # area and hole status are invariant under rotation
  return(p)
}

rotate <- function(X, ...) {
  UseMethod("rotate")
}

rotate.owin <- function(X, angle=pi/2, ..., rescue=TRUE, centre=NULL) {
  verifyclass(X, "owin")
  if(!is.null(centre)) {
    ## rotation about designated centre
    X <- shift(X, origin=centre)
    negorig <- getlastshift(X)
  } else negorig <- NULL
  switch(X$type,
         rectangle={
           # convert rectangle to polygon
           P <- owin(X$xrange, X$yrange, poly=
                     list(x=X$xrange[c(1,2,2,1)],
                          y=X$yrange[c(1,1,2,2)]),
                     unitname=unitname(X))
           # call polygonal case
           Y <- rotate.owin(P, angle, rescue=rescue)
         },
         polygonal={
           # First rotate the polygonal boundaries
           bdry <- lapply(X$bdry, rotxypolygon, angle=angle)
           # wrap up
           Y <- owin(poly=bdry, check=FALSE, unitname=unitname(X))
           if(rescue)
             Y <- rescue.rectangle(Y)
         },
         mask={
           newframe <- boundingbox(rotxy(corners(X), angle))
           Y <- if(length(list(...)) > 0) as.mask(newframe, ...) else 
                   as.mask(newframe, eps=with(X, min(xstep, ystep)))
           pixelxy <- rasterxy.mask(Y)
           xybefore <- rotxy(pixelxy, -angle)
           Y$m[] <- with(xybefore, inside.owin(x, y, X))
           Y <- intersect.owin(Y, boundingbox(Y))
           if(rescue)
             Y <- rescue.rectangle(Y)
           unitname(Y) <- unitname(X)
         },
         stop("Unrecognised window type")
         )
  if(!is.null(negorig))
    Y <- shift(Y, -negorig)
  return(Y)
}

rotate.ppp <- function(X, angle=pi/2, ..., centre=NULL) {
  verifyclass(X, "ppp")
  if(!is.null(centre)) {
    X <- shift(X, origin=centre)
    negorigin <- getlastshift(X)
  } else negorigin <- NULL
  r <- rotxy(X, angle)
  w <- rotate.owin(X$window, angle, ...)
  Y <- ppp(r$x, r$y, window=w, marks=marks(X, dfok=TRUE), check=FALSE)
  if(!is.null(negorigin))
    Y <- shift(Y, -negorigin)
  return(Y)
}

rotate.im <- function(X, angle=pi/2, ..., centre=NULL) {
  if(!is.null(centre)) {
    X <- shift(X, origin=centre)
    negorigin <- getlastshift(X)
  } else negorigin <- NULL
  co <- cos(angle)
  si <- sin(angle)
  m <- matrix(c(co,si,-si,co), nrow=2, ncol=2)
  Y <- affine(X, mat=m)
  if(!is.null(negorigin))
    Y <- shift(Y, -negorigin)
  return(Y)
}

