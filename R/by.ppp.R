#
#   by.ppp.R
#
#  $Revision: 1.6 $  $Date: 2015/10/21 09:06:57 $
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
  z <- as.solist(z, demote=TRUE)
  return(z)
}
