#
# plot.mppm.R
#
#   $Revision: 1.3 $  $Date: 2015/09/27 02:28:17 $
#
#

plot.mppm <- function(x, ..., trend=TRUE, cif=FALSE, se=FALSE,
                      how=c("image", "contour", "persp")) {
  xname <- deparse(substitute(x))
  how <- match.arg(how)
  subs <- subfits(x)
  arglist <- resolve.defaults(list(x=subs, how=how),
                              list(...),
                              list(main=xname))
  if(trend) 
    do.call("plot", c(arglist, list(trend=TRUE, cif=FALSE, se=FALSE)))
  if(cif) 
    do.call("plot", c(arglist, list(trend=FALSE, cif=TRUE, se=FALSE)))
  if(se) 
    do.call("plot", c(arglist, list(trend=FALSE, cif=FALSE, se=TRUE)))
  invisible(NULL)
}

