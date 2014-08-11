#
# plot.mppm.R
#
#   $Revision: 1.1 $  $Date: 2007/03/26 01:31:40 $
#
#

plot.mppm <- function(x, ..., trend=TRUE, cif=FALSE, how="image") {
  xname <- deparse(substitute(x))
  if(length(how) > 1)
    stop(paste("Multiple plotting styles cannot be selected;",
               "argument", dQuote("how"), "must have length 1"))
  if(!missing(trend) && missing(cif))
    cif <- !trend
  else if(missing(trend) && !missing(cif))
    trend <- !cif
  else if(trend + cif != 1)
    stop(paste("Exactly one of", dQuote("trend"), "and", dQuote("cif"),
               "should be TRUE"))
  subs <- subfits(x)
  arglist <- resolve.defaults(list(x=subs,how=how, trend=trend, cif=cif),
                              list(...),
                              list(main=xname))
  do.call("plot", arglist)
}

