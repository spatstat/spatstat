#
#   ssf.R
#
#  spatially sampled functions
#
#  $Revision: 1.17 $  $Date: 2017/01/26 00:55:22 $
#

ssf <- function(loc, val) {
  stopifnot(is.ppp(loc))
  if(is.function(val))
    val <- val(loc$x, loc$y)
  if(is.data.frame(val))
    val <- as.matrix(val)
  if(!is.matrix(val))
    val <- matrix(val, ncol=1, dimnames=list(NULL, "value"))
  if(nrow(val) != npoints(loc))
    stop("Incompatible lengths")
  result <- loc %mark% val
  class(result) <- c("ssf", class(result))
  attr(result, "ok") <- complete.cases(val)
  return(result)
}

print.ssf <- function(x, ..., brief=FALSE) {
  if(brief) {
    cat(paste("Spatial function sampled at", npoints(x), "locations\n"))
  } else {
    cat("Spatially sampled function\n")
    cat("Locations:\n\t")
    print(unmark(x))
  }
  val <- marks(x)
  if(!is.matrix(val)) {
    d <- 1
    warning("Internal format error: val is not a matrix")
  } else d <- ncol(val) 
  if(!brief) {
    type <- if(d == 1) "Scalar" else paste(d, "-vector", sep="")
    cat(paste(type, "valued function\n"))
  }
  if(d > 1 && !is.null(nama <- colnames(val)))
    cat(paste("Component names:", commasep(sQuote(nama)), "\n"))
  return(invisible(NULL))
}

image.ssf <- function(x, ...) {
  do.call("plot", resolve.defaults(list(x, how="smoothed"), list(...)))
}

as.im.ssf <- function(X, ...) nnmark(X, ...)

as.function.ssf <- function(x, ...) {
  X <- x
  mX <- marks(X)
  switch(markformat(X),
         vector = {
           g <- function(x, y=NULL) {
             Y <- xy.coords(x,y)[c("x","y")]
             J <- nncross(Y, X, what="which")
             result <- mX[J]
             return(unname(result))
           }
         },
         dataframe = {
           g <- function(x, y=NULL) {
             Y <- xy.coords(x,y)[c("x","y")]
             J <- nncross(Y, X, what="which")
             result <-  mX[J,,drop=FALSE]
             row.names(result) <- NULL
             return(result)
           }
         },
         stop("Marks must be a vector or data.frame"))
  h <- funxy(g, Frame(X))
  return(h)
}

plot.ssf <- function(x, ..., how=c("smoothed", "nearest", "points"),
                     style = c("image", "contour", "imagecontour"),
                     sigma=NULL, contourargs=list()) {
  xname <- short.deparse(substitute(x))
  how <- match.arg(how)
  style <- match.arg(style)
  otherargs <- list(...)
  # convert to images
  y <- switch(how,
              points = as.ppp(x),
              nearest = nnmark(x), 
              smoothed = Smooth(x, sigma=sigma)
              )
  # points plot
  if(how == "points") {
    out <- do.call("plot",
                   resolve.defaults(list(y), otherargs,
                                    list(main=xname)))
    if(is.null(out)) return(invisible(NULL))
    return(out)
  }
  # image plot
  switch(style,
         image = {
           out <- do.call("plot",
                          resolve.defaults(list(y), otherargs,
                                           list(main=xname)))
         },
         contour = {
           do.call("plot",
                   resolve.defaults(list(as.owin(x)),
                                    otherargs, list(main=xname)))
           do.call("contour",
                   resolve.defaults(list(y, add=TRUE), contourargs))
           out <- NULL
         },
         imagecontour = {
           out <- do.call("plot",
                          resolve.defaults(list(y), otherargs,
                                           list(main=xname)))
           do.call("contour",
                   resolve.defaults(list(y, add=TRUE), contourargs))
         })
  return(invisible(out))
}

contour.ssf <- function(x, ..., main, sigma=NULL) {
  if(missing(main))
    main <- short.deparse(substitute(x))
  y <- Smooth(x, sigma=sigma)
  contour(y, ..., main=main)
  return(invisible(NULL))
}

Smooth.ssf <- function(X, ...) {
  stopifnot(inherits(X, "ssf"))
  ok  <- attr(X, "ok")
  Y   <- as.ppp(X)[ok]
  argh <- list(...)
  isnul <- as.logical(unlist(lapply(argh, is.null)))
  nonnularg <- argh[!isnul]
  sigma0 <- if(any(c("sigma", "varcov") %in% names(nonnularg)))
            NULL else 1.4 * max(nndist(X))
  Z <- do.call("Smooth.ppp",
               resolve.defaults(list(X = Y),
                                list(...),
                                list(sigma=sigma0),
                                .MatchNull=FALSE))
                                # don't take NULL for an answer!
  return(Z)
}

"[.ssf" <-
  function(x, i, j, ..., drop) {
  loc <- unmark(x)
  val <- marks(x)
  ok  <- attr(x, "ok")
  #
  if(!missing(j)) 
    val <- val[, j, drop=FALSE]
  if(!missing(i)) {
    # use [.ppp to identify which points are retained
    locn <- loc %mark% seq_len(npoints(loc))
    loci <- locn[i]
    loc  <- unmark(loci)
    id   <- marks(loci)
    # extract
    val  <- val[id, , drop=FALSE]
    ok   <- ok[id]
  }
  out <- loc %mark% val
  class(out) <- c("ssf", class(out))
  attr(out, "ok") <- ok
  return(out)    
}

as.ppp.ssf <- function(X, ...) {
  class(X) <- "ppp"
  attr(X, "ok") <- NULL
  return(X)
}

marks.ssf <-  function(x, ...) {
  val <- x$marks
  if(is.null(dim(val))) val <- matrix(val, ncol=1)
  if(is.data.frame(val)) val <- as.matrix(val)
  return(val)
}

"marks<-.ssf" <- function(x, ..., value) {
  ssf(unmark(x), value)
}

unmark.ssf <- function(X) { unmark(as.ppp(X)) }

with.ssf <- function(data, ...) {
  loc <- as.ppp(data)
  val <- marks(data)
  newval <- with(as.data.frame(val), ...)
  if(length(newval) == npoints(loc) ||
     (is.matrix(newval) && nrow(newval) == npoints(loc)))
    return(ssf(loc, newval))
  return(newval)
}

apply.ssf <- function(X, ...) {
  loc <- as.ppp(X)
  val <- marks(X)
  newval <- apply(val, ...)
  if(length(newval) == npoints(loc) ||
     (is.matrix(newval) && nrow(newval) == npoints(loc)))
    return(ssf(loc, newval))
  return(newval)
}

range.ssf <- function(x, ...) range(marks(x), ...)
min.ssf <- function(x, ...) min(marks(x), ...)
max.ssf <- function(x, ...) max(marks(x), ...)

integral.ssf <- function(f, domain=NULL, ..., weights=attr(f, "weights")) {
  if(!is.null(weights)) {
    check.nvector(weights, npoints(f), oneok=TRUE)
    if(length(weights) == 1) weights <- rep(weights, npoints(f))
  }
  if(!is.null(domain)) {
    ok <- inside.owin(f, w=domain)
    f <- f[ok,]
    if(!is.null(weights)) weights <- weights[ok]
  }
  y <- marks(f)
  if(is.null(weights)) {
    z <- if(!is.matrix(y)) mean(y, na.rm=TRUE) else colMeans(y, na.rm=TRUE)
    a <- area(Window(f))
  } else {
    z <- if(!is.matrix(y)) weighted.mean(y, w=weights, na.rm=TRUE) else 
         apply(y, 2, weighted.mean, w=weights, na.rm=TRUE)
    a <- sum(weights)
  }
  return(z * a)
}
