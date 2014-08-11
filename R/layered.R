#
# layered.R
#
# Simple mechanism for layered plotting
#
#  $Revision: 1.24 $  $Date: 2014/01/15 06:57:15 $
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
    np <- length(plotargs)
    if(np == 1) plotargs <- rep(plotargs, n) else if(np != n)
      stop("plotargs should have one component for each element of the list")
  } else {
    plotargs <- rep.int(list(list()), n)
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

plot.layered <- function(x, ..., which=NULL, plotargs=NULL,
                         add=FALSE, show.all=!add, main=NULL) {
  if(is.null(main))
    main <- short.deparse(substitute(x))
  n <- length(x)
  if(!is.null(plotargs)) {
    np <- length(plotargs)
    if(!(is.list(plotargs) && all(unlist(lapply(plotargs, is.list)))))
      stop("plotargs should be a list of lists")
  }
  ## select layers
  if(!is.null(which)) {
    x <- x[which]
    nw <- length(x)
    if(!is.null(plotargs)) {
      if(np == n) plotargs <- plotargs[which] else
      if(np == 1) plotargs <- rep(plotargs, nw) else
      if(np != nw) 
        stop("plotargs should have one component for each layer to be plotted")
    }
    n <- nw
  } else if(!is.null(plotargs)) {
    if(np == 1) plotargs <- rep(plotargs, n) else
    if(np != n) stop("plotargs should have one component for each layer")
  }
  ## remove null layers
  if(any(isnul <- unlist(lapply(x, is.null)))) {
    x <- x[!isnul]
    if(!is.null(plotargs))
      plotargs <- plotargs[!isnul]
    n <- length(x)
  }
  ## anything to plot?
  if(n == 0)
    return(invisible(NULL))
  ## Merge plotting arguments
  xplotargs <- layerplotargs(x)
  if(is.null(plotargs)) {
    plotargs <- xplotargs
  } else if(length(xplotargs) > 0) {
    for(i in 1:n)
      plotargs[[i]] <- resolve.defaults(plotargs[[i]], xplotargs[[i]])
  }
  ## Determine bounding box 
  bb <- NULL
  boxes <- lapply(x,
                  function(z) { try(as.rectangle(z), silent=TRUE) })
  if(!any(unlist(lapply(boxes, inherits, what="try-error")))) 
    bb <- do.call("bounding.box", boxes)
  ## Start plotting
  started <- inked <- FALSE
  if(add) {
    started <- TRUE
    if(show.all && !is.null(bb)) 
      fakemaintitle(bb, main, ...)
  } else if(!is.null(bb)) {
    ## initialise new plot using bounding box
    bb <- do.call("bounding.box", boxes)
    plot(bb, type="n", main=if(show.all) main else "")
    started <- TRUE
  }
  # plot the layers
  out <- list()
  nama <- names(x)
  for(i in seq_along(x)) {
    xi <- x[[i]]
    if(length(xi) == 0) {
      # null layer - no plotting
      out[[i]] <- NULL
    } else {
      ## plot layer i on top of previous layers if any.
      ## Show all graphic elements of the first component only;
      ## but do not display the names of any components.
      show.name.i <- resolve.1.default(list(show.all=FALSE),
                                       list(...), 
                                       plotargs[[i]])
      dflt <- list(main=if(show.name.i) nama[i] else "",
                   show.all=show.all && !inked)
      ## 
      out[[i]] <- do.call("plot",
                          resolve.defaults(list(x=xi,
                                                add=started),
                                           list(...),
                                           plotargs[[i]],
                                           dflt))
      started <- titled <- TRUE
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
  if(!(is.list(value) && all(unlist(lapply(value, is.list)))))
    stop("Replacement value should be a list of lists")
  n <- length(L)
  if(length(value) == 1) value <- unname(rep(value, n)) else 
  if(length(value) != n) stop("Replacement value is wrong length")
  if(is.null(names(value))) names(value) <- names(L) else
  if(!identical(names(value), names(L)))
    stop("Mismatch in names of list elements")
  attr(L, "plotargs") <- value
  return(L)
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
  
shift.layered <- function(X, vec=c(0,0), ...) {
  if(length(list(...)) > 0) {
    if(!missing(vec)) 
      warning("Argument vec ignored; overridden by other arguments")
    ## ensure the same shift is applied to all layers
    s <- shift(X[[1]], ...)
    vec <- getlastshift(s)
  }
  Y <- applytolayers(X, shift, vec=vec)
  attr(Y, "lastshift") <- vec
  return(Y)
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
  if(length(W) == 0) {
    if(fatal) stop("Layered object has no window data")
    return(NULL)
  }
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

