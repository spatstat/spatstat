##
## solist.R
##
## Methods for class `solist' (spatial object list)
##
##      and related classes 'anylist', 'ppplist', 'imlist'
##
## plot.solist is defined in plot.solist.R
##
## $Revision: 1.5 $ $Date: 2014/10/18 01:43:45 $

anylist <- function(...) {
  x <- list(...)
  class(x) <- c("anylist", "listof", class(x))
  return(x)
}

print.anylist <- function (x, ...) {
  nn <- names(x)
  ll <- length(x)
  if (length(nn) != ll) 
    nn <- paste("Component", seq.int(ll))
  for (i in seq_len(ll)) {
    splat(paste0(nn[i], ":"))
    print(x[[i]], ...)
    cat("\n")
  }
  return(invisible(x))
}

as.anylist <- function(x) {
  if(inherits(x, "anylist")) return(x)
  if(!is.list(x))
    x <- list(x)
  class(x) <- c("anylist", "listof", class(x))
  return(x)
}
  
"[.anylist" <- function(x, i, ...) {
  cl <- oldClass(x)
  ## invoke list method
  y <- NextMethod("[")
  if(length(y) == 0) return(list())
  class(y) <- cl
  return(y)
}

"[<-.anylist" <- function(x, i, value) {
  as.anylist(NextMethod("[<-"))
}

summary.anylist <- function(object, ...) {
  as.anylist(lapply(object, summary, ...))
}

## .................... solist .............................

solist <- local({

  checkspatial <-
    function(z) !inherits(try(Frame(z), silent=TRUE), "try-error")

  solist <- function(..., check=TRUE, promote=TRUE) {
    stuff <- list(...)
    if(check && !all(unlist(lapply(stuff, checkspatial))))
      stop("Some arguments of solist() are not 2D spatial objects")
    class(stuff) <- c("solist", "anylist", "listof", class(stuff))
    if(promote) {
      if(all(unlist(lapply(stuff, is.ppp)))) {
        class(stuff) <- c("ppplist", class(stuff))
      } else if(all(unlist(lapply(stuff, is.im)))) {
        class(stuff) <- c("imlist", class(stuff))
      }
    }
    return(stuff)
  }

  solist
})

as.solist <- function(x, ...) {
  if(inherits(x, "solist") && length(list(...)) == 0)
    return(x)
  if(!is.list(x))
    x <- list(x)
  return(do.call(solist, append(x, list(...))))
}

print.solist <- function (x, ...) {
  what <- if(inherits(x, "ppplist")) "point patterns" else
          if(inherits(x, "imlist")) "pixel images" else "spatial objects"
  splat(paste("List of", what))
  parbreak()
  NextMethod("print")
}


"[.solist" <- function(x, i, ...) {
  cl <- oldClass(x)
  ## invoke list method
  y <- NextMethod("[")
  if(length(y) == 0) return(list())
  class(y) <- cl
  return(y)
}
  
"[<-.solist" <- function(x, i, value) {
  ## invoke list method
  y <- NextMethod("[<-")
  ## check again
  return(do.call(solist, y))
}
  
summary.solist <- function(object, ...) {
  x <- lapply(object, summary, ...)
  attr(x, "otype") <-
    if(inherits(object, "ppplist")) "ppp" else
    if(inherits(object, "imlist")) "im" else ""
  class(x) <- c("summary.solist", "anylist")
  x
}

print.summary.solist <- function(x, ...) {
  what <- switch(attr(x, "otype"),
                 ppp="point patterns",
                 im="pixel images",
                 "spatial objects")
  splat("Summary of", length(x), what)
  parbreak()
  NextMethod("print")
}

as.layered.solist <- function(X) {
  layered(LayerList=X)
}

