#
# layered.R
#
# Simple mechanism for layered plotting
#
#  $Revision: 1.15 $  $Date: 2013/02/25 08:45:04 $
#

layered <- function(..., plotargs=NULL, LayerList=NULL) {
  argh <- list(...)
  if(length(argh) > 0 && !is.null(LayerList))
    stop("LayerList is incompatible with other arguments")
  out <- if(!is.null(LayerList)) LayerList else argh
  n <- length(out)
  if(sum(nzchar(names(out))) != n)
    names(out) <- paste("Layer", seq_len(n))
  if(!is.null(plotargs)) {
    if(!is.list(plotargs) || !all(unlist(lapply(plotargs, is.list))))
      stop("plotargs should be a list of lists")
    if(length(plotargs) != length(out))
      stop("plotargs should have one component for each element of the list")
  } else {
    plotargs <- rep(list(list()), length(out))
  }
  names(plotargs) <- names(out)
  attr(out, "plotargs") <- plotargs
  class(out) <- c("layered", class(out))
  return(out)
}

print.layered <- function(x, ...) {
  cat("Layered object\n")
  for(i in seq_along(x)) {
    cat(paste("\n", names(x)[i], ":\n", sep=""))
    print(x[[i]])
  }
  pl <- layerplotargs(x)
  hasplot <- (unlist(lapply(pl, length)) > 0)
  if(any(hasplot)) 
    cat(paste("\nIncludes plot arguments for",
              commasep(names(pl)[hasplot]), "\n"))
  invisible(NULL)
}

plot.layered <- function(x, ..., which=NULL, plotargs=NULL) {
  xname <- short.deparse(substitute(x))
  xp <- if(is.null(which)) x else x[which]
  if(length(xp) == 0)
    return(invisible(NULL))
  # validate plotting arguments
  if(is.null(plotargs)) {
    plotargs <- attr(x, "plotargs")
    if(!is.null(plotargs) && !is.null(which)) plotargs <- plotargs[which]
  } else {
    if(!is.list(plotargs) || !all(unlist(lapply(plotargs, is.list))))
      stop("plotargs should be a list of lists")
    if(length(plotargs) != length(xp))
      stop("plotargs should have one component for each layer to be plotted")
  }
  # determine plot frame 
  add <- resolve.1.default("add", list(...), list(add=FALSE))
  if(add) {
    started <- TRUE
  } else {
    # new plot
    # determine bounding frame
    isnul <- unlist(lapply(x, is.null))
    boxes <- lapply(x[!isnul],
                    function(z) { try(as.rectangle(z), silent=TRUE) })
    if(any(unlist(lapply(boxes, inherits, what="try-error")))) {
       # failed to determine bounding frame
      started <- FALSE
    } else {
      bb <- do.call("bounding.box", boxes)
      plot(bb, type="n", main=xname)
      started <- TRUE
    }
  }
  # plot the layers
  out <- list()
  for(i in seq_along(xp)) {
    xpi <- xp[[i]]
    if(length(xpi) == 0) {
      # null layer - no plotting
      out[[i]] <- NULL
    } else {
      # plot layer i on top of previous layers
      iargs <- if(!started) list(main=xname) else list(add=TRUE)
      out[[i]] <- do.call("plot",
                          resolve.defaults(list(x=xpi),
                                           list(...),
                                           plotargs[[i]],
                                           iargs))
      started <- TRUE
    }
  }
  return(invisible(out))
}

"[.layered" <- function(x, i, j, drop=FALSE, ...) {
  if(missing(i) && missing(j))
    return(x)
  p <- attr(x, "plotargs")
  x <- unclass(x)
  nx <- length(x)
  if(!missing(i) && !is.null(i)) {
    x <- x[i]
    p <- p[i]
    nx <- length(x)
  }
  isnul <- (unlist(lapply(x, length)) == 0)
  if(!missing(j) && !is.null(j))
    x[!isnul] <- lapply(x[!isnul], "[", i=j)
  if(drop && nx == 1)
    return(x[[1]])
  y <- layered(LayerList=x, plotargs=p)
  return(y)
}

layerplotargs <- function(L) {
  stopifnot(inherits(L, "layered"))
  attr(L, "plotargs")
}

"layerplotargs<-" <- function(L, value) {
  stopifnot(inherits(L, "layered"))
  if(length(value) != length(L))
    stop("Replacement value is wrong length")
  if(!identical(names(value), names(L)))
    stop("Mismatch in names of list elements")
  attr(L, "plotargs") <- value
  L
}

applytolayers <- function(L, FUN, ...) {
  # Apply FUN to each **non-null** layer,
  # preserving the plot arguments
  pla <- layerplotargs(L)
  ok <- !unlist(lapply(L, is.null))
  L[ok] <- lapply(L[ok], FUN, ...)
  Z <- layered(LayerList=L, plotargs=pla)
  return(Z)
}
  
shift.layered <- function(X, ...) {
  applytolayers(X, shift, ...)
}

affine.layered <- function(X, ...) {
  applytolayers(X, affine, ...)
}

rotate.layered <- function(X, ...) {
  applytolayers(X, rotate, ...)
}

reflect.layered <- function(X) {
  applytolayers(X, reflect)
}

flipxy.layered <- function(X) {
  applytolayers(X, flipxy)
}

scalardilate.layered <- function(X, ...) {
  applytolayers(X, scalardilate, ...)
}
  
rescale.layered <- function(X, s) {
  if(!missing(s)) applytolayers(X, rescale, s=s) else applytolayers(X, rescale)
}


as.owin.layered <- function(W, ..., fatal=TRUE) {
  if(length(W) == 0) {
    if(fatal) stop("Layered object is empty: no window data")
    return(NULL)
  }
  # remove null layers
  isnul <- unlist(lapply(W, is.null))
  W <- W[!isnul]
  Wlist <- lapply(unname(W), as.owin, ..., fatal=fatal)
  Wlist <- lapply(Wlist, rescue.rectangle)
  Z <- Wlist[[1]]
  if(length(Wlist) > 1) {
    same <- unlist(lapply(Wlist[-1], identical, y=Z))
    if(!all(same))
      Z <- do.call("union.owin", Wlist)
  }
  return(Z)
}
