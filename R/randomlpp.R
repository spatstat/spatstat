#
#  random.R
#
#  Random point pattern generators for a linear network
#
#  $Revision: 1.4 $   $Date: 2015/02/17 06:15:44 $
#

rpoislpp <- function(lambda, L, ..., nsim=1) {
  if(missing(L) || is.null(L)) {
    if(!inherits(lambda, c("linim", "linfun")))
      stop("L is missing", call.=FALSE)
    L <- as.linnet(lambda)
  } else verifyclass(L, "linnet")
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
