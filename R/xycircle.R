#'
#'   xycircle.R
#'
#'   Low-level utilities for circle geometry
#'
#'  $Revision: 1.8 $   $Date: 2017/01/21 10:57:04 $
#'

xysegXcircle <- function(xcentres, ycentres, radii, x0, y0, x1, y1,
                         check=TRUE) {
  #' 'Cross' version
  #' Find intersections between circles and segments
  #' for all combinations of centres, radii and segments.
  #'
  #'   xcentres, ycentres: numeric vectors of coordinates of centres
  #'   radii:              numeric vector of radii
  #'   x0, y0, x1, y1:     numeric vectors of segment endpoint coordinates
  #'
  if(check)
    stopifnot(is.numeric(xcentres),
              is.numeric(ycentres),
              length(xcentres) == length(ycentres),
              is.numeric(radii),
              is.numeric(x0),
              is.numeric(y0),
              is.numeric(x1),
              is.numeric(y1),
              length(x0) == length(y0),
              length(x1) == length(y1),
              length(x0) == length(x1))
  storage.mode(xcentres) <- storage.mode(ycentres) <- "double"
  storage.mode(x0) <- storage.mode(y0) <- "double"
  storage.mode(x1) <- storage.mode(y1) <- "double"
  storage.mode(radii) <- "double"
  z <- .Call("circXseg",
             XC = xcentres,
             YC = ycentres,
             R  = radii,
             X0 = x0,
             Y0 = y0,
             X1 = x1,
             Y1 = y1)
  result <- as.data.frame(z)
  #' indices i, j, k specify provenance of each intersection point
  #' i = circle, j = segment, k = radius
  names(result) <- c("x", "y", "i", "j", "k", "sinalpha")
  return(result)
}

xysegMcircle <- function(xcentres, ycentres, radmat, x0, y0, x1, y1,
                         check=TRUE) {
  #' 'Matrix' version
  #' Find intersections between circles and segments
  #' where radii are given in a matrix with rows corresponding to centres.
  #'
  #'   xcentres, ycentres: numeric vectors of coordinates of centres
  #'   radmat:             matrix of radii (rows correspond to centres)
  #'   x0, y0, x1, y1:     numeric vectors of segment endpoint coordinates
  #'
  if(check)
    stopifnot(is.numeric(xcentres),
              is.numeric(ycentres),
              length(xcentres) == length(ycentres),
              is.numeric(radmat),
              is.matrix(radmat),
              nrow(radmat) == length(xcentres),
              is.numeric(x0),
              is.numeric(y0),
              is.numeric(x1),
              is.numeric(y1),
              length(x0) == length(y0),
              length(x1) == length(y1),
              length(x0) == length(x1))
  storage.mode(xcentres) <- storage.mode(ycentres) <- "double"
  storage.mode(x0) <- storage.mode(y0) <- "double"
  storage.mode(x1) <- storage.mode(y1) <- "double"
  storage.mode(radmat) <- "double"
  z <- .Call("circMseg",
             XC = xcentres,
             YC = ycentres,
             R  = radmat,
             X0 = x0,
             Y0 = y0,
             X1 = x1,
             Y1 = y1)
  result <- as.data.frame(z)
  #' indices i, j, k specify provenance of each intersection point
  #' i = circle, j = segment, k = radius
  names(result) <- c("x", "y", "i", "j", "k", "sinalpha")
  return(result)
}

xysegPcircle <- function(xc, yc, rc, x0, y0, x1, y1,
                         check=TRUE) {
  #' 'Parallel' version
  #' Find intersections between circles and segments
  #' for circles with centres (xc, yc) and radii (rc) corresponding.
  #'
  #'   xc, y:   numeric vectors of coordinates of centres
  #'   rc:      numeric vector of radii (corresponding to xc, yc)
  #'   x0, y0, x1, y1:    numeric vectors of segment endpoint coordinates
  #'
  if(check)
    stopifnot(is.numeric(xc),
              is.numeric(yc),
              length(xc) == length(yc),
              is.numeric(rc),
              length(rc) == length(xc),
              is.numeric(x0),
              is.numeric(y0),
              is.numeric(x1),
              is.numeric(y1),
              length(x0) == length(y0),
              length(x1) == length(y1),
              length(x0) == length(x1))
  storage.mode(xc) <- storage.mode(yc) <- "double"
  storage.mode(x0) <- storage.mode(y0) <- "double"
  storage.mode(x1) <- storage.mode(y1) <- "double"
  storage.mode(rc) <- "double"
  z <- .Call("circXseg",
             XC = xc,
             YC = yc,
             RC = rc,
             X0 = x0,
             Y0 = y0,
             X1 = x1,
             Y1 = y1)
  result <- as.data.frame(z)
  #' indices i, j specify provenance of each intersection point
  #' i = circle, j = segment
  names(result) <- c("x", "y", "i", "j", "sinalpha")
  return(result)
}

