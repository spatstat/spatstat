#
# plot.mppm.R
#
#   $Revision: 1.6 $  $Date: 2020/12/19 05:25:06 $
#
#

plot.mppm <- function(x, ..., trend=TRUE, cif=FALSE, se=FALSE,
                      how=c("image", "contour", "persp")) {
  xname <- deparse(substitute(x))
  how <- match.arg(how)
  subs <- subfits(x)
  dont.complain.about(subs)
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

