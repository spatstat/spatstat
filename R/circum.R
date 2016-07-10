#'
#'    circum.R
#'
#'  circumradius, circumcentre and circumdisc
#'
#'  $Revision: 1.2 $ $Date: 2016/07/10 03:22:06 $
#'

circumradius <- function(x, ...) {
  UseMethod("circumradius")
}

circumcentre <- function(x, ...) {
  UseMethod("circumcentre")
}

circumcircle <- function(x, ...) {
  UseMethod("circumcircle")
}

#' owin

circumradius.owin <- function(x, ...) {
  sqrt(min(fardist(x, ..., squared=TRUE)))
}

circumcentre.owin <- function(x, ...) {
  z <- where.min(fardist(x, ..., squared=TRUE))
  Window(z) <- x
  return(z)
}

circumcircle.owin <- function(x, ...) {
  d2 <- fardist(x, ..., squared=TRUE)
  z <- where.min(d2)
  r <- sqrt(min(d2))
  w <- disc(centre=z, radius=r) 
  return(w)
}

#' ppp

circumradius.ppp <- function(x, ...) {
  circumradius(convexhull(x), ...)
}

circumcentre.ppp <- function(x, ...) {
  z <- circumcentre(convexhull(x), ...)
  Window(z) <- Window(x)
  return(z)
}

circumcircle.ppp <- function(x, ...) {
  circumcircle(convexhull(x), ...)
}

