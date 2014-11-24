#
# listof.R
#
# Methods for class `listof'
#
# plot.listof is defined in plot.splitppp.R
#

"[<-.listof" <- function(x, i, value) {
  # invoke list method
  class(x) <- "list"
  x[i] <- value
  # then make it a 'listof' object too
  class(x) <- c("listof", class(x))
  x
}
  
summary.listof <- function(object, ...) {
  x <- lapply(object, summary, ...)
  class(x) <- "summary.listof"
  x
}

print.summary.listof <- function(x, ...) {
  class(x) <- "listof"
  print(x)
  invisible(NULL)
}

listof <- function(...) {
#  warn.once("listof",
#            "The class listof will be Deprecated",
#            "in future versions of spatstat.",
#            "Use anylist or solist")
  stuff <- list(...)
  class(stuff) <- c("listof", class(stuff))
  return(stuff)
}

as.listof <- function(x) {
  if(!is.list(x))
    x <- list(x)
  if(!inherits(x, "listof"))
    class(x) <- c("listof", class(x))
#  warn.once("listof",
#            "The class listof will be Deprecated",
#            "in future versions of spatstat.",
#            "Use anylist or solist")
  return(x)
}

contour.listof <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  do.call("plot.listof",
          resolve.defaults(list(x=x, plotcommand="contour"),
                           list(...),
                           list(main=xname)))
}

image.listof <- local({

  image.listof <- function(x, ..., equal.ribbon = FALSE, ribmar=NULL) {
    xname <- short.deparse(substitute(x))
    if(equal.ribbon && !all(unlist(lapply(x, is.im)))) {
      warning("equal.ribbon is only implemented for objects of class 'im'")
      equal.ribbon <- FALSE
    }
    if(equal.ribbon) imagecommon(x, ..., xname=xname, ribmar=ribmar) else 
      do.call("plot.listof",
              resolve.defaults(list(x=x, plotcommand="image"),
                               list(...),
                               list(main=xname)))
  }

  imagecommon <- function(x, ...,
                          xname,
                          zlim=NULL,
                          ribbon=TRUE,
                          ribside=c("right", "left", "bottom", "top"),
                          ribsep=NULL, ribwid=0.5, ribn=1024,
                          ribscale=NULL, ribargs=list(),
                          ribmar=NULL,
                          mar.panel=c(2,1,1,2)) {
    if(missing(xname))
      xname <- short.deparse(substitute(x))
    ribside <- match.arg(ribside)
    stopifnot(is.list(ribargs))
    if(!is.null(ribsep))
      warning("Argument ribsep is not yet implemented for image arrays")
    # determine range of values
    if(is.null(zlim))
      zlim <- range(unlist(lapply(x, range)))
    # determine common colour map
    imcolmap <- plot.im(x[[1]], do.plot=FALSE, zlim=zlim, ..., ribn=ribn)
    # plot ribbon?
    if(!ribbon) {
      ribadorn <- list()
    } else {
      # determine plot arguments for colour ribbon
      vertical <- (ribside %in% c("right", "left"))
      scaleinfo <- if(!is.null(ribscale)) list(labelmap=ribscale) else list()
      sidecode <- match(ribside, c("bottom", "left", "top", "right"))
      ribstuff <- c(list(x=imcolmap, main="", vertical=vertical),
                    ribargs,
                    scaleinfo,
                    list(side=sidecode))
      if(is.null(ribmar)) {
          ribmar <- mar.panel/2
          newmar <- c(2, 0) ## 2 for title, 0 for opposite side
          switch(ribside,
                 left   = { ribmar[c(2,4)] <- newmar },
                 right  = { ribmar[c(4,2)] <- newmar },
                 bottom = { ribmar[c(1,3)] <- newmar },
                 top    = { ribmar[c(3,1)] <- newmar })
      }
      # function executed to plot colour ribbon
      do.ribbon <- function() {
        opa <- par(mar=ribmar)
        do.call("plot", ribstuff)
        par(opa)
      }
      # encoded as 'adorn' argument
      ribadorn <- list(adorn=do.ribbon, adorn.size=ribwid)
      names(ribadorn)[1] <- paste("adorn", ribside, sep=".")
    }
    #
    do.call("plot.listof",
            resolve.defaults(list(x=x, plotcommand="image"),
                             list(...),
                             list(main=xname),
                             list(col=imcolmap, zlim=zlim, ribbon=FALSE),
                             ribadorn))
  }

  image.listof
})

as.layered.listof <- function(X) {
  layered(LayerList=X)
}
