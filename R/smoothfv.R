#
#  smoothfv.R
#
#   $Revision: 1.4 $   $Date: 2010/04/16 12:47:43 $
#
  
smooth.fv <- function(x, which="*", ...) {
  stopifnot(is.character(which))
  if(length(which) == 1 && which %in% .Spatstat.FvAbbrev) {
    if(which == ".x")
      stop("Cannot smooth the function argument")
    which <- fvnames(x, which)
  }
  if(any(nbg <- !(which %in% names(x)))) 
    stop(paste("Unrecognised column",
               ngettext(sum(nbg), "name", "names"),
               commasep(sQuote(which[nbg])), 
               "in argument", sQuote("which")))
  xx <- x[[fvnames(x, ".x")]]
  for(ynam in which) {
    yy <- x[[ynam]]
    ok <- is.finite(yy)
    ss <- smooth.spline(xx[ok], yy[ok], ...)
    yhat <- predict(ss, xx)$y
    x[[ynam]] <- yhat
  }
  return(x)
}
