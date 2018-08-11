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

as.layered.listof <- function(X) {
  layered(LayerList=X)
}
