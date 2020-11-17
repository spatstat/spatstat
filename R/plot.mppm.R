#
# plot.mppm.R
#
#   $Revision: 1.5 $  $Date: 2020/11/17 03:47:24 $
#
#

plot.mppm <- function(x, ..., trend=TRUE, cif=FALSE, se=FALSE,
                      how=c("image", "contour", "persp")) {
  xname <- deparse(substitute(x))
  how <- match.arg(how)
  subs <- subfits(x)
  arglist <- resolve.defaults(list(x=quote(subs), how=how),
                              list(...),
                              list(main=xname))
  if(trend) 
    do.call(plot, c(arglist, list(trend=TRUE, cif=FALSE, se=FALSE)))
  if(cif) 
    do.call(plot, c(arglist, list(trend=FALSE, cif=TRUE, se=FALSE)))
  if(se) 
    do.call(plot, c(arglist, list(trend=FALSE, cif=FALSE, se=TRUE)))
  invisible(NULL)
}

