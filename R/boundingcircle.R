#'
#'    boundingcircle.R
#'
#'  bounding circle and its centre
#'
#'  $Revision: 1.5 $ $Date: 2016/07/16 03:07:02 $
#'

circumradius <- function(x, ...) {
  .Deprecated("boundingradius")
  UseMethod("boundingradius")
}
circumradius.owin <- function(x, ...) {
  .Deprecated("boundingradius.owin")
  boundingradius.owin(x, ...)
}
circumradius.ppp <- function(x, ...) {
  .Deprecated("boundingradius.ppp")
  boundingradius.ppp(x, ...)
}

boundingradius <- function(x, ...) {
  UseMethod("boundingradius")
}

boundingcentre <- function(x, ...) {
  UseMethod("boundingcentre")
}

boundingcircle <- function(x, ...) {
  UseMethod("boundingcircle")
}

#' owin

boundingradius.owin <- function(x, ...) {
  sqrt(min(fardist(x, ..., squared=TRUE)))
}

boundingcentre.owin <- function(x, ...) {
  z <- where.min(fardist(x, ..., squared=TRUE))
  Window(z) <- x
  return(z)
}

boundingcircle.owin <- function(x, ...) {
  d2 <- fardist(x, ..., squared=TRUE)
  z <- where.min(d2)
  r <- sqrt(min(d2))
  w <- disc(centre=z, radius=r) 
  return(w)
}

#' ppp

boundingradius.ppp <- function(x, ...) {
  boundingradius(convexhull(x), ...)
}

boundingcentre.ppp <- function(x, ...) {
  z <- boundingcentre(convexhull(x), ...)
  Window(z) <- Window(x)
  return(z)
}

boundingcircle.ppp <- function(x, ...) {
  boundingcircle(convexhull(x), ...)
}

