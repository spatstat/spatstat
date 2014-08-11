#
#  smoothfv.R
#
#   $Revision: 1.8 $   $Date: 2013/08/05 07:08:51 $
#
  
smooth.fv <- function(x, which="*", ..., 
                      method=c("smooth.spline", "loess"),
                      xinterval=NULL) {
  stopifnot(is.character(which))
  method <- match.arg(method)
  if(!is.null(xinterval))
    check.range(xinterval) 
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
  # process each column of function values
  for(ynam in which) {
    yy <- x[[ynam]]
    ok <- is.finite(yy)
    if(!is.null(xinterval))
      ok <- ok & inside.range(xx, xinterval)
    switch(method,
           smooth.spline = {
             ss <- smooth.spline(xx[ok], yy[ok], ...)
             yhat <- predict(ss, xx[ok])$y
           },
           loess = {
             df <- data.frame(x=xx[ok], y=yy[ok])
             lo <- loess(y ~ x, df, ...)
             yhat <- predict(lo, df[,"x", drop=FALSE])
           })
    yy[ok] <- yhat
    x[[ynam]] <- yy
  }
  return(x)
}
