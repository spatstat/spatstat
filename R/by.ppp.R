#
#   by.ppp.R
#
#  $Revision: 1.5 $  $Date: 2011/05/18 01:29:48 $
#

by.ppp <- function(data, INDICES=marks(data), FUN, ...) {
  if(missing(INDICES))
    INDICES <- marks(data, dfok=FALSE)
  if(missing(FUN))
    stop("FUN is missing")
  y <- split(data, INDICES)
  z <- list()
  for(i in seq_along(y))
    z[[i]] <- FUN(y[[i]], ...)
  names(z) <- names(y)
  z <- as.listof(z)
  return(z)
}
