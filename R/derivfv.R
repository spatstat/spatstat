#
# derivfv.R
#
# differentiation for fv objects
#
#  $Revision: 1.4 $ $Date: 2013/08/05 07:08:31 $
#

deriv.fv <- function(expr, which="*", ...,
                     method=c("spline", "numeric"),
                     kinks=NULL) {
  f <- expr
  method <- match.arg(method)
  # select columns
  if(length(which) == 1 && which %in% .Spatstat.FvAbbrev) {
    if(which == ".x")
      stop("Cannot smooth the function argument")
    which <- fvnames(f, which)
  }
  if(any(nbg <- !(which %in% names(f)))) 
    stop(paste("Unrecognised column",
               ngettext(sum(nbg), "name", "names"),
               commasep(sQuote(which[nbg])), 
               "in argument", sQuote("which")))
  relevant <- names(f) %in% which
  # get 
  rname <- fvnames(f, ".x")
  df <- as.data.frame(f)
  rpos <- which(colnames(df) == rname)
  rvals <- df[,rpos]
  yvals <- df[,relevant,drop=FALSE]
  nr <- length(rvals)
  # cut x axis into intervals?
  if(is.null(kinks)) {
    cutx <- factor(rep(1, nr))
  } else {
    rr <- range(rvals)
    breaks <- sort(unique(kinks))
    if(breaks[1] > rr[1]) breaks <- c(rr[1], breaks)
    if(max(breaks) < rr[2]) breaks <- c(breaks, rr[2])
    cutx <- cut(rvals, breaks=breaks, include.lowest=TRUE)
  }
  # process
  for(segment in levels(cutx)) {
    ii <- (cutx == segment)
    yy <- yvals[ii, , drop=FALSE]
    switch(method,
           numeric = {
             dydx <- apply(yy, 2, diff)/diff(rvals[ii])
             nd <- nrow(dydx)
             dydx <- rbind(dydx, dydx[nd, ])
         },
         spline = {
           dydx <- apply(yy, 2,
                         function(y, r, ...) {
                           ss <- smooth.spline(r, y, ...)
                           predict(ss, r, deriv=1)$y
                         },
                         r=rvals[ii], ...)
         })
    df[ii, relevant] <- dydx
  }
  desc <- attr(expr, "desc")
  desc[relevant] <- paste("derivative of", desc[relevant])
  result <- fv(df,
               argu=rname,
               ylab=NULL,
               valu=fvnames(expr, ".y"),
               fmla=attr(expr, "fmla"),
               alim=attr(expr, "alim"),
               desc=desc,
               unitname=unitname(expr),
               fname=paste("D", attr(expr, "fname"), sep=""))
  fvnames(result, ".") <- fvnames(f, ".")
  return(result)
}
