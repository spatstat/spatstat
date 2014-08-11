#
#	rotate.S
#
#	$Revision: 1.18 $	$Date: 2012/10/10 01:20:23 $
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

rotate.owin <- function(X, angle=pi/2, ..., rescue=TRUE) {
  verifyclass(X, "owin")
  switch(X$type,
         rectangle={
           # convert rectangle to polygon
           P <- owin(X$xrange, X$yrange, poly=
                     list(x=X$xrange[c(1,2,2,1)],
                          y=X$yrange[c(1,1,2,2)]),
                     unitname=unitname(X))
           # call polygonal case
           return(rotate.owin(P, angle, rescue=rescue))
         },
         polygonal={
           # First rotate the polygonal boundaries
           bdry <- lapply(X$bdry, rotxypolygon, angle=angle)
           # wrap up
           W <- owin(poly=bdry, check=FALSE, unitname=unitname(X))
           if(rescue)
             W <- rescue.rectangle(W)
           return(W)
         },
         mask={
           newframe <- bounding.box.xy(rotxy(corners(X), angle))
           W <- if(length(list(...)) > 0) as.mask(newframe, ...) else 
                   as.mask(newframe, eps=with(X, min(xstep, ystep)))
           pixelxy <- raster.xy(W)
           xybefore <- rotxy(pixelxy, -angle)
           W$m[] <- with(xybefore, inside.owin(x, y, X))
           W <- intersect.owin(W, bounding.box(W))
           if(rescue)
             W <- rescue.rectangle(W)
           unitname(W) <- unitname(X)
           return(W)
         },
         stop("Unrecognised window type")
         )
}

rotate.ppp <- function(X, angle=pi/2, ...) {
  verifyclass(X, "ppp")
  r <- rotxy(X, angle)
  w <- rotate.owin(X$window, angle, ...)
  return(ppp(r$x, r$y, window=w, marks=marks(X, dfok=TRUE), check=FALSE))
}

rotate.im <- function(X, angle=pi/2, ...) {
  co <- cos(angle)
  si <- sin(angle)
  m <- matrix(c(co,si,-si,co), nrow=2, ncol=2)
  affine(X, mat=m)
}

