#
# layered.R
#
# Simple mechanism for layered plotting
#
#  $Revision: 1.34 $  $Date: 2015/10/21 09:06:57 $
#

layered <- function(..., plotargs=NULL, LayerList=NULL) {
  argh <- list(...)
  if(length(argh) > 0 && !is.null(LayerList))
    stop("LayerList is incompatible with other arguments")
  out <- if(!is.null(LayerList)) LayerList else argh
  n <- length(out)
  if(sum(nzchar(names(out))) != n)
    names(out) <- paste("Layer", seq_len(n))
  if(is.null(plotargs)) {
    plotargs <- rep.int(list(list()), n)
  } else {
    if(!is.list(plotargs))
      stop("plotargs should be a list of lists")
    if(!all(unlist(lapply(plotargs, is.list))))
      plotargs <- list(plotargs)
    np <- length(plotargs)
    if(np == 1) plotargs <- rep(plotargs, n) else if(np != n)
      stop("plotargs should have one component for each element of the list")
  }
  names(plotargs) <- names(out)
  attr(out, "plotargs") <- plotargs
  class(out) <- c("layered", class(out))
  return(out)
}

print.layered <- function(x, ...) {
  splat("Layered object")
  if(length(x) == 0) splat("(no entries)")
  for(i in seq_along(x)) {
    cat(paste("\n", names(x)[i], ":\n", sep=""))
    print(x[[i]])
  }
  pl <- layerplotargs(x)
  hasplot <- (unlist(lapply(pl, length)) > 0)
  if(any(hasplot)) 
    splat("Includes plot arguments for", commasep(names(pl)[hasplot]))
  invisible(NULL)
}

plot.layered <- function(x, ..., which=NULL, plotargs=NULL,
                         add=FALSE, show.all=!add, main=NULL,
                         do.plot=TRUE) {
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
  a <- plotEachLayer(x, ..., plotargs=plotargs, add=add,
                     show.all=show.all, do.plot=FALSE)
  if(!do.plot)
    return(a)
  bb <- as.rectangle(as.owin(a))
  ## Start plotting
  if(!add && !is.null(bb)) {
    ## initialise new plot using bounding box
    pt <- prepareTitle(main)
    plot(bb, type="n", main=pt$blank)
    add <- TRUE
  }
  # plot the layers
  out <- plotEachLayer(x, ..., main=main,
                       plotargs=plotargs, add=add,
                       show.all=show.all, do.plot=TRUE)
  return(invisible(out))
}

plotEachLayer <- function(x, ..., main,
                          plotargs, add, show.all, do.plot=TRUE) {
  main.given <- !missing(main)
  ## do.plot=TRUE    =>   plot the layers 
  ## do.plot=FALSE   =>   determine bounding boxes
  out <- boxes <- list()
  nama <- names(x)
  firstlayer <- TRUE
  for(i in seq_along(x)) {
    xi <- x[[i]]
    if(length(xi) == 0) {
      # null layer - no plotting
      out[[i]] <- boxes[[i]] <- NULL
    } else {
      ## plot layer i on top of previous layers if any.
      ## By default,
      ##    - show all graphic elements of the first component only;
      ##    - show title 'firstmain' on first component;
      ##    - do not show any component names.
      add.i <- add || !firstlayer
      if(main.given) {
        main.i <- if(firstlayer) main else ""
      } else {
        show.all.i <- resolve.1.default(list(show.all=FALSE),
                                         list(...), 
                                         plotargs[[i]])
        main.i <- if(show.all.i) nama[i] else ""
      }
      dflt <- list(main=main.i,
                   show.all=show.all && firstlayer)
      pla.i <- plotargs[[i]]
      defaultplot <- !(".plot" %in% names(pla.i))
      ## plot layer i, or just determine bounding box
      if(defaultplot &&
         inherits(xi, c("ppp", "psp", "owin",
                        "lpp", "linnet", 
                        "im", "msr", "layered"))) {
        ## plot method for 'xi' has argument 'do.plot'.
        out[[i]] <- outi <- do.call("plot",
                                    resolve.defaults(list(x=xi,
                                                          add=add.i,
                                                          do.plot=do.plot),
                                                     list(...),
                                                     pla.i,
                                                     dflt))
        boxes[[i]] <- as.rectangle(as.owin(outi))
      } else {
        ## plot method for 'xi' does not have argument 'do.plot'
        if(do.plot) {
          if(defaultplot) {
            plotfun <- "plot"
          } else {
            plotfun <- pla.i[[".plot"]]
            pla.i <- pla.i[names(pla.i) != ".plot"]
          }
          out[[i]] <- outi <- do.call(plotfun,
                                      resolve.defaults(list(x=xi,
                                                            add=add.i),
                                                       list(...),
                                                       pla.i,
                                                       dflt))
        }
        ## convert layer i to box
        boxi <- try(as.rectangle(xi), silent=TRUE)
        boxes[[i]] <- if(!inherits(boxi, "try-error")) boxi else NULL
      }
      firstlayer <- FALSE
    }
  }
  ## one box to bound them all
  if(!all(unlist(lapply(boxes, is.null))))
    attr(out, "bbox") <- do.call(boundingbox, boxes)
  return(out)
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

"[[<-.layered" <- function(x, i, value) {
  x[i] <- if(!is.null(value)) list(value) else NULL
  return(x)
}

"[<-.layered" <- function(x, i, value) {
  p <- layerplotargs(x)
  ## invoke list method
  y <- x
  class(y) <- "list"
  y[i] <- value
  # make it a 'layered' object too
  class(y) <- c("layered", class(y))
  # update names and plotargs
  if(any(blank <- !nzchar(names(y)))) {
    names(y)[blank] <- paste("Layer", which(blank))
    pnew <- rep(list(list()), length(y))
    names(pnew) <- names(y)
    m <- match(names(y), names(x))
    mok <- !is.na(m)
    pnew[mok] <- p[m[mok]]
    layerplotargs(y) <- pnew
  } else layerplotargs(y) <- layerplotargs(x)[names(y)]
  return(y)
}

layerplotargs <- function(L) {
  stopifnot(inherits(L, "layered"))
  attr(L, "plotargs")
}

"layerplotargs<-" <- function(L, value) {
  if(!inherits(L, "layered"))
    L <- layered(L)
  if(!is.list(value))
    stop("Replacement value should be a list, or a list-of-lists")
  n <- length(L)
  if(!all(unlist(lapply(value, is.list)))) 
    value <- unname(rep(list(value), n))
  if(length(value) != n) {
    if(length(value) == 1) value <- unname(rep(value, n)) else
    stop("Replacement value is wrong length")
  }
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
  if(length(L) > 0) {
    ok <- !unlist(lapply(L, is.null))
    L[ok] <- lapply(L[ok], FUN, ...)
  }
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

rotate.layered <- function(X, ..., centre=NULL) {
  if(!is.null(centre)) {
    X <- shift(X, origin=centre)
    negorigin <- getlastshift(X)
  } else negorigin <- NULL
  Y <- applytolayers(X, rotate, ...)
  if(!is.null(negorigin))
    Y <- shift(Y, -negorigin)
  return(Y)
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
  
rescale.layered <- function(X, s, unitname) {
  if(missing(s)) s <- NULL
  if(missing(unitname)) unitname <- NULL
  applytolayers(X, rescale, s=s, unitname=unitname) 
}


as.owin.layered <- local({

  as.owin.layered <- function(W, ..., fatal=TRUE) {
    if(length(W) == 0) {
      if(fatal) stop("Layered object is empty: no window data")
      return(NULL)
    }
    ## remove null layers
    isnul <- unlist(lapply(W, is.null))
    W <- W[!isnul]
    if(length(W) == 0) {
      if(fatal) stop("Layered object has no window data")
      return(NULL)
    }
    Wlist <- lapply(unname(W), as.owin, ..., fatal=fatal)
    Wlist <- lapply(Wlist, rescue.rectangle)
    Wlist <- lapply(Wlist, puffbox)
    Z <- Wlist[[1]]
    if(length(Wlist) > 1) {
      same <- unlist(lapply(Wlist[-1], identical, y=Z))
      if(!all(same))
        Z <- do.call("union.owin", Wlist)
    }
    return(Z)
  }

  puffbox <- function(W) {
    ## union.owin will delete boxes that have width zero or height zero
    ## so 'puff' them out slightly
    ss <- sidelengths(Frame(W))
    if(ss[1] == 0) W$xrange <- W$xrange + 1e-6 * c(-1,1) * ss[2]
    if(ss[2] == 0) W$yrange <- W$yrange + 1e-6 * c(-1,1) * ss[1]
    return(W)
  }
  
  as.owin.layered
})


domain.layered <- Window.layered <- function(X, ...) { as.owin(X) }

as.layered <- function(X) {
  UseMethod("as.layered")
}

as.layered.default <- function(X) {
  layered(X)
}

as.layered.ppp <- function(X) {
  if(!is.marked(X)) return(layered(X))
  if(is.multitype(X)) return(layered(LayerList=split(X)))
  mX <- marks(X)
  if(!is.null(d <- dim(mX)) && d[2] > 1) {
    mx <- as.data.frame(marks(X))
    Y <- lapply(mx, setmarks, x=X)
    return(layered(LayerList=Y))
  }
  return(layered(X))
}


  
