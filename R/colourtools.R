#
#  colourtools.R
#
#   $Revision: 1.17 $   $Date: 2015/07/08 10:23:36 $
#


rgb2hex <- function(v, maxColorValue=255) {
  stopifnot(is.numeric(v))
  if(!is.matrix(v))
    v <- matrix(v, nrow=1)
  if(ncol(v) %in% c(3, 4)) {
    out <- rgb(v, maxColorValue=maxColorValue)
    return(out)
  } 
  stop("v should be a vector of length 3 or 4, or a matrix with 3 or 4 columns")
}

rgb2hsva <- function(red, green=NULL, blue=NULL, alpha=NULL,
                     maxColorValue=255) {
  if(is.null(green) && is.null(blue) && is.null(alpha)) {
    ## red should be a 3-row matrix of RGB values
    ## or a 4-row matrix of RGBA values 
    if(!is.matrix(red))
      red <- matrix(red, ncol=1)
    ## check for an alpha channel
    if(nrow(red) == 4) {
      alpha <- red[4,]
      red <- red[-4, , drop=FALSE]
    }
  }
  y <- rgb2hsv(red, green, blue, maxColorValue=maxColorValue)
  if(!is.null(alpha))
    y <- rbind(y, alpha=alpha/maxColorValue)
  return(y)
}
 
col2hex <- function(x) {
  # convert to RGBA
  y <- col2rgb(x, alpha=TRUE)
  # remove alpha channel if all colours are opaque
  if(all(y["alpha", ] == 255))
    y <- y[1:3, , drop=FALSE]
  # convert to hex 
  z <- rgb2hex(t(y))
  return(z)
}

paletteindex <- function(x) {
  x <- col2hex(x)
  p <- col2hex(palette())
  m <- match(x, p)
  return(m)
}

is.colour <- function(x) {
  if(length(x) == 0) return(FALSE)
  cx <- try(col2rgb(x), silent=TRUE)
  bad <- inherits(cx, "try-error")
  return(!bad)
}

samecolour <- function(x, y) { col2hex(x) == col2hex(y) }

complementarycolour <- function(x) {
  if(is.null(x)) return(NULL)
  if(inherits(x, "colourmap")) {
    colouroutputs(x) <- complementarycolour(colouroutputs(x))
    return(x)
  }
  # convert to RGBA
  y <- col2rgb(x, alpha=TRUE)
  # complement of R, G, B
  y[1:3, ] <- 255 - y[1:3, ]
  # convert to colours
  z <- rgb2hex(t(y))
  return(z)
}

is.grey <- function(x) {
  if(inherits(x, "colourmap")) x <- colouroutputs(x)
  if(is.function(x)) return(NA)
  y <- rgb2hsva(col2rgb(x, alpha=TRUE))
  sat <- y["s", ]
  alp <- y["alpha", ]
  return(sat == 0 & alp == 1)
}

to.opaque <- function(x) {
  if(all(!is.na(paletteindex(x))))
    return(x) # preserve palette colours
  rgb(t(col2rgb(x)), maxColorValue=255)
}

to.transparent <- function(x, fraction) {
  if(all(fraction == 1))
    return(to.opaque(x))
  rgb(t(col2rgb(x))/255, alpha=fraction, maxColorValue=1)
}
  
to.grey <- function(x, weights=c(0.299, 0.587, 0.114), transparent=FALSE) {
  if(is.null(x)) return(NULL)
  if(inherits(x, "colourmap")) {
    colouroutputs(x) <- to.grey(colouroutputs(x),
                                weights=weights, transparent=transparent)
    return(x)
  }
  if(is.function(x)) {
    f <- x
    g <- function(...) to.grey(f(...), weights=weights, transparent=transparent)
    return(g)
  }
  ## preserve palette indices, if only using black/grey
  if(all(!is.na(paletteindex(x))) && all(is.grey(x)))
    return(x)
  if(!transparent) {
    y <- col2rgb(x)
    g <- (weights %*% y)/(255 * sum(weights))
    z <- grey(g)
  } else {
    yy <- col2rgb(x, alpha=TRUE)
    y <- yy[1:3, , drop=FALSE]
    g <- (weights %*% y)/(255 * sum(weights))
    z <- grey(g, alpha=y[4,])
  }
  return(z)
}

is.col.argname <- function(x) {
  return(nzchar(x) & ((x == "col") | (substr(x, 1, 4) == "col.")))
}

col.args.to.grey <- function(x, ...) {
  if(any(hit <- is.col.argname(names(x))))
    x[hit] <- lapply(x[hit], to.grey, ...)
  return(x)
}

# versions of rgb() and hsv() that work with NA values

rgbNA <- function(red, green, blue, alpha=NULL, maxColorValue=1) {
  df <- if(is.null(alpha)) data.frame(red=red, green=green, blue=blue) else
        data.frame(red=red, green=green, blue=blue, alpha=alpha)
  result <- rep(NA_character_, nrow(df))
  ok <- complete.cases(df)
  result[ok] <- if(is.null(alpha)) {
    with(df, rgb(red[ok], green[ok], blue[ok],
                 maxColorValue=maxColorValue))
  } else {
    with(df, rgb(red[ok], green[ok], blue[ok], alpha[ok],
                 maxColorValue=maxColorValue))
  }
  return(result)
}

hsvNA <- function(h, s, v, alpha=NULL) {
  df <- if(is.null(alpha)) data.frame(h=h, s=s, v=v) else
                           data.frame(h=h, s=s, v=v, alpha=alpha)
  result <- rep(NA_character_, nrow(df))
  ok <- complete.cases(df)
  result[ok] <- if(is.null(alpha)) {
    with(df, hsv(h[ok], s[ok], v[ok]))
  } else {  
    with(df, hsv(h[ok], s[ok], v[ok], alpha[ok]))
  }
  return(result)
}

