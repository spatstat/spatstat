#
#  random.R
#
#  Random point pattern generators for a linear network
#
#  $Revision: 1.3 $   $Date: 2012/10/20 06:56:01 $
#

rpoislpp <- function(lambda, L, ...) {
  verifyclass(L, "linnet")
  X <- datagen.rpoisppOnLines(lambda, as.psp(L), ...)
  Y <- lpp(X, L)
  return(Y)
}

runiflpp <- function(n, L) {
  verifyclass(L, "linnet")
  X <- datagen.runifpointOnLines(n, as.psp(L))
  Y <- lpp(X, L)
  return(Y)
}
