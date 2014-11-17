#
#  random.R
#
#  Random point pattern generators for a linear network
#
#  $Revision: 1.3 $   $Date: 2012/10/20 06:56:01 $
#

rpoislpp <- function(lambda, L, ..., nsim=1) {
  verifyclass(L, "linnet")
  result <- vector(mode="list", length=nsim)
  for(i in 1:nsim) {
    X <- datagen.rpoisppOnLines(lambda, as.psp(L), ...)
    Y <- lpp(X, L)
    if(nsim == 1) return(Y)
    result[[i]] <- Y
  }
  Y <- as.solist(Y)
  names(Y) <- paste("Simulation", 1:nsim)
  return(Y)
}

runiflpp <- function(n, L, nsim=1) {
  verifyclass(L, "linnet")
  result <- vector(mode="list", length=nsim)
  for(i in 1:nsim) {
    X <- datagen.runifpointOnLines(n, as.psp(L))
    Y <- lpp(X, L)
    if(nsim == 1) return(Y)
    result[[i]] <- Y
  }
  Y <- as.solist(Y)
  names(Y) <- paste("Simulation", 1:nsim)
  return(Y)
}
