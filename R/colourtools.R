#
#  colourtools.R
#
#   $Revision: 1.10 $   $Date: 2015/01/19 12:00:33 $
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

complementarycolour <- function(x) {
  if(is.null(x)) return(NULL)
  if(inherits(x, "colourmap")) {
    colouroutputs(x) <- complementarycolour(colouroutputs(x))
    return(x)
  }
  y <- apply(255 - col2rgb(x), 2, rgb2hex)
  return(y)
}

is.grey <- function(x) {
  if(inherits(x, "colourmap")) x <- colouroutputs(x)
  if(is.function(x)) return(NA)
  sat <- rgb2hsv(col2rgb(x))["s", ]
  return(sat == 0)
}
  
to.grey <- function(x, weights=c(0.299, 0.587, 0.114)) {
  if(is.null(x)) return(NULL)
  if(inherits(x, "colourmap")) {
    colouroutputs(x) <- to.grey(colouroutputs(x), weights=weights)
    return(x)
  }
  if(is.function(x)) {
    f <- x
    g <- function(...) to.grey(f(...))
    return(g)
  }
  ## preserve palette indices, if only using black/grey
  if(all(!is.na(paletteindex(x))) && all(is.grey(x)))
    return(x)
  y <- col2rgb(x)
  z <- (weights %*% y)/(255 * sum(weights))
  return(grey(z))
}

# versions of rgb() and hsv() that work with NA values

rgbNA <- function(red, green, blue, ...) {
  with(data.frame(red=red, green=green, blue=blue), {
    ok <- !(is.na(red) | is.na(green) | is.na(blue))
    values <- rgb(red[ok], green[ok], blue[ok], ...)
    result <- character(length(red))
    result[ok] <- values
    result[!ok] <- NA
    return(result)
  })
}

hsvNA <- function(h, s, v, ...) {
  with(data.frame(h=h, s=s, v=v), {
    ok <- !(is.na(h) | is.na(s) | is.na(v))
    values <- hsv(h[ok], s[ok], v[ok], ...)
    result <- character(length(h))
    result[ok] <- values
    result[!ok] <- NA
    return(result)
  })
}

