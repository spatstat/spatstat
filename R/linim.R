#
# linim.R
#
#  $Revision: 1.7 $   $Date: 2013/01/15 04:57:55 $
#
#  Image/function on a linear network
#

linim <- function(L, Z, ..., df=NULL) {
  L <- as.linnet(L)
  stopifnot(is.im(Z))
  if(is.null(df)) {
    # compute the data frame of mapping information
    xx <- rasterx.im(Z)
    yy <- rastery.im(Z)
    mm <- !is.na(Z$v)
    xx <- as.vector(xx[mm])
    yy <- as.vector(yy[mm])
    pixelcentres <- ppp(xx, yy, window=as.rectangle(Z), check=FALSE)
    pixdf <- data.frame(xc=xx, yc=yy)
    # project pixel centres onto lines
    p2s <- project2segment(pixelcentres, as.psp(L))
    projloc <- as.data.frame(p2s$Xproj)
    projmap <- as.data.frame(p2s[c("mapXY", "tp")])
    # extract values
    values <- Z[pixelcentres]
    # bundle
    df <- cbind(pixdf, projloc, projmap, data.frame(values=values))
  } else {
    stopifnot(is.data.frame(df))
    neednames <- c("xc", "yc", "x", "y", "mapXY", "tp", "values")
    ok <- neednames %in% names(df)
    if(any(!ok)) {
      nn <- sum(!ok)
      stop(paste(ngettext(nn, "A column", "Columns"),
                 "named", commasep(sQuote(neednames[!ok])),
                 ngettext(nn, "is", "are"),
                 "missing from argument", sQuote("df")))
    }
  }
  out <- Z
  attr(out, "L") <- L
  attr(out, "df") <- df
  class(out) <- c("linim", class(out))
  return(out)
}

print.linim <- function(x, ...) {
  cat("Image on linear network\n")
  print(attr(x, "L"))
  NextMethod("print")
}

plot.linim <- function(x, ..., style=c("colour", "width"), scale, adjust=1) {
  xname <- short.deparse(substitute(x))
  style <- match.arg(style)
  # colour style: plot as pixel image
  if(style == "colour")
    return(do.call("plot.im",
                   resolve.defaults(list(x),
                                    list(...),
                                    list(main=xname))))
  # width style
  L <- attr(x, "L")
  df <- attr(x, "df")
  Llines <- as.psp(L)
  # initialise plot
  W <- as.owin(L)
  do.call.matched("plot.owin",
                  resolve.defaults(list(x=W, type="n"),
                                   list(...), list(main=xname)),
                  extrargs="type")
  # rescale values to a plottable range
  vr <- range(df$values)
  vr[1] <- min(0, vr[1])
  if(missing(scale)) {
    maxsize <- mean(distmap(Llines))/2
    scale <- maxsize/diff(vr)
  } 
  df$values <- adjust * scale * (df$values - vr[1])
  # split data by segment
  mapXY <- factor(df$mapXY, levels=seq_len(Llines$n))
  dfmap <- split(df, mapXY, drop=TRUE)
  # sort each segment's data by position along segment
  dfmap <- lapply(dfmap, function(z) { z[fave.order(z$tp), ] })
  # plot each segment's data
  Lends <- Llines$ends
  Lperp <- angles.psp(Llines) + pi/2
  Lfrom <- L$from
  Lto   <- L$to
  Lvert <- L$vertices
  for(i in seq(length(dfmap))) {
    z <- dfmap[[i]]
    segid <- unique(z$mapXY)[1]
    xx <- z$x
    yy <- z$y
    vv <- z$values
    # add endpoints of segment
    leftend <- Lvert[Lfrom[segid]]
    rightend <- Lvert[Lto[segid]]
    xx <- c(leftend$x, xx, rightend$x)
    yy <- c(leftend$y, yy, rightend$y)
    vv <- c(vv[1],     vv, vv[length(vv)])
    # create polygon
    xx <- c(xx, rev(xx))
    yy <- c(yy, rev(yy))
    vv <- c(vv, -rev(vv))/2
    ang <- Lperp[segid]
    xx <- xx + cos(ang) * vv
    yy <- yy + sin(ang) * vv
    do.call.matched("polygon",
                    resolve.defaults(list(x=xx, y=yy),
                                     list(...),
                                     list(border=NA, col=1)))
  }
  return(invisible(NULL))
}

as.im.linim <- function(X, ...) { as.im(X$Z, ...) }

as.linim <- function(X, ...) {
  UseMethod("as.linim")
}

as.linim.default <- function(X, L, ...) {
  stopifnot(inherits(L, "linnet"))
  Y <- as.im(X, W=as.rectangle(as.owin(L)), ...)
  Z <- as.im(as.mask.psp(as.psp(L), as.owin(Y)))
  Y <- eval.im(Z * Y)
  out <- linim(L, Y)
  return(out)
}

as.linim.linim <- function(X, ...) {
  if(length(list(...)) == 0)
    return(X)
  Y <- as.linim.default(X, as.linnet(X), ...)
  return(Y)
}

    
