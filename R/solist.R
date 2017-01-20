##
## solist.R
##
## Methods for class `solist' (spatial object list)
##
##      and related classes 'anylist', 'ppplist', 'imlist'
##
## plot.solist is defined in plot.solist.R
##
## $Revision: 1.14 $ $Date: 2017/01/20 10:06:37 $

anylist <- function(...) {
  x <- list(...)
  class(x) <- c("anylist", "listof", class(x))
  return(x)
}

print.anylist <- function (x, ...) {
  ll <- length(x)
  if(ll == 0) {
    splat("(Zero length list)")
    return(invisible(NULL))
  }
  nn <- names(x)
  if (length(nn) != ll) 
    nn <- paste("Component", seq.int(ll))
  spaceok <- waxlyrical('space')
  for (i in seq_len(ll)) {
    splat(paste0(nn[i], ":"))
    print(x[[i]], ...)
    if(spaceok && i < ll) cat("\n")
  }
  return(invisible(NULL))
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

pool.anylist <- function(x, ...) {
  do.call(pool, append(x, list(...)))
}

## .................... solist .............................

is.sob <- local({
  ## test whether x is a spatial object suitable for solist
  sobjectclasses <- c("ppp", "psp", "im", "owin", 
                      "quad", "tess", "msr",
                      "quadratcount", "quadrattest", 
                      "layered",
                      "funxy", "distfun", "nnfun", 
                      "lpp", "linnet", "linfun",      
                      "influence.ppm", "leverage.ppm")
  # Note 'linim' inherits 'im'
  #      'dfbetas.ppm' inherits 'msr'

  is.sob <- function(x) { inherits(x, what=sobjectclasses) }
  is.sob
})
  
solist <- function(..., check=TRUE, promote=TRUE, demote=FALSE) {
  stuff <- list(...)
  if((check || demote) && !all(sapply(stuff, is.sob))) {
    if(demote)
      return(as.anylist(stuff))
    stop("Some arguments of solist() are not 2D spatial objects")
  }
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

as.solist <- function(x, ...) {
  if(inherits(x, "solist") && length(list(...)) == 0)
    return(x)
  if(!is.list(x) || is.sob(x))
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
  if(!missing(i) && is.owin(i)) {
    ## spatial subset
    y <- lapply(unclass(x), "[", i=i)
  } else {
    ## invoke list method
    y <- NextMethod("[")
  }
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

# --------------- counterparts of 'lapply' --------------------

anylapply <- function(X, FUN, ...) {
  v <- lapply(X, FUN, ...)
  return(as.anylist(v))
}

solapply <- function(X, FUN, ..., check=TRUE, promote=TRUE, demote=FALSE) {
  v <- lapply(X, FUN, ...)
  u <- as.solist(v, check=check, promote=promote, demote=demote)
  return(u)
}

density.ppplist <- function(x, ..., se=FALSE) {
  y <- lapply(x, density, ..., se=se)
  if(!se) return(as.solist(y, demote=TRUE))
  y.est <- lapply(y, getElement, name="estimate")
  y.se  <- lapply(y, getElement, name="SE")
  z <- list(estimate = as.solist(y.est, demote=TRUE),
            SE       = as.solist(y.se,  demote=TRUE))
  return(z)
}

expandSpecialLists <- function(x, special="solist") {
  ## x is a list which may include entries which are lists, of class 'special'
  ## unlist these entries only
  hit <- sapply(x, inherits, what=special) 
  if(!any(hit)) return(x)
  # wrap each *non*-special entry in list()
  x[!hit] <- lapply(x[!hit], list)
  # now strip one layer of list() from all entries
  return(unlist(x, recursive=FALSE))
}
