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
  stuff <- list(...)
  class(stuff) <- c("listof", class(stuff))
  return(stuff)
}

as.listof <- function(x) {
  if(!is.list(x))
    x <- list(x)
  if(!inherits(x, "listof"))
    class(x) <- c("listof", class(x))
  return(x)
}

contour.listof <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  do.call("plot.listof",
          resolve.defaults(list(x=x, plotcommand="contour"),
                           list(...),
                           list(main=xname)))
}

image.listof <- function(x, ..., equal.ribbon = FALSE) {
  xname <- short.deparse(substitute(x))
  if(equal.ribbon) {
    dotargs <- list(...)
    # all panels will be coded using the same colour map
    # determine range of values
    zlim <- dotargs$zlim
    # compute colour map 
    if(!is.null(zlim)) {
      imcolmap <- plot.im(x[[1]], preponly=TRUE, ...)
    } else {
      zlim <- range(unlist(lapply(x, range)))
      imcolmap <- plot.im(x[[1]], preponly=TRUE, zlim=zlim, ...)
    }
    # assemble arguments for plot.listof
    rightribbon <- function() {
      opa <- par(mar=c(1,0,1,2))
      plot(imcolmap, vertical=TRUE, main="")
      par(opa)
    }
    
    zz <- list(col=imcolmap, zlim=zlim,
               ribbon=FALSE,
               adorn.right=rightribbon)
  } else zz <- list()
  do.call("plot.listof",
          resolve.defaults(list(x=x, plotcommand="image"),
                           zz,
                           list(...),
                           list(main=xname)))
}
