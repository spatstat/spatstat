#
#  colourtools.R
#
#   $Revision: 1.6 $   $Date: 2014/01/29 03:59:19 $
#


rgb2hex <- function(v) {
  stopifnot(is.numeric(v))
  if(is.matrix(v)) {
    stopifnot(ncol(v) == 3)
  } else {
    if(length(v) != 3)
      stop("v should be a vector of length 3 or a matrix with 3 columns")
    v <- matrix(v, ncol=3)
  } 
  out <- rgb(v[,1], v[,2], v[,3], maxColorValue=255)
  return(out)
}

col2hex <- function(x) { apply(col2rgb(x), 2, rgb2hex) }

paletteindex <- function(x) {
  x <- col2hex(x)
  p <- col2hex(palette())
  m <- match(x, p)
  return(m)
}

samecolour <- function(x, y) { col2hex(x) == col2hex(y) }

to.grey <- function(x, weights=c(1,1,1)) {
  if(is.null(x)) return(NULL)
  if(inherits(x, "colourmap")) {
    col <- summary(x)[["outputs"]]
    return(tweak.colourmap(x, to.grey(col, weights)))
  }
  y <- col2rgb(x)
  z <- (weights %*% y)/(255 * sum(weights))
  return(grey(z))
}

# versions of rgb() and hsv() that work with NA values

rgbNA <- function(red, green, blue, ...) {
  ok <- !(is.na(red) | is.na(green) | is.na(blue))
  values <- rgb(red[ok], green[ok], blue[ok], ...)
  result <- character(length(red))
  result[ok] <- values
  result[!ok] <- NA
  return(result)
}

hsvNA <- function(h, s, v, ...) {
  ok <- !(is.na(h) | is.na(s) | is.na(v))
  values <- hsv(h[ok], s[ok], v[ok], ...)
  result <- character(length(h))
  result[ok] <- values
  result[!ok] <- NA
  return(result)
}

