#'
#'   xycircle.R
#'
#'   Low-level utilities for circle geometry
#'
#'  $Revision: 1.4 $   $Date: 2017/01/20 11:44:29 $
#'

xysegXcircle <- function(xcentres, ycentres, radii, x0, y0, x1, y1,
                         check=TRUE) {
  #' find intersections between circles and segments
  #' (all combinations of centres, radii and segments are considered)
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
  storage.mode(r) <- "double"
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
