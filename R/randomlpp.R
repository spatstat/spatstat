#
#  random.R
#
#  Random point pattern generators for a linear network
#
#  $Revision: 1.5 $   $Date: 2016/10/31 00:49:56 $
#

rpoislpp <- function(lambda, L, ..., nsim=1, drop=TRUE) {
  if(missing(L) || is.null(L)) {
    if(!inherits(lambda, c("linim", "linfun")))
      stop("L is missing", call.=FALSE)
    L <- as.linnet(lambda)
  } else verifyclass(L, "linnet")
  result <- vector(mode="list", length=nsim)
  S <- as.psp(L)
  bugout <- (nsim == 1) && drop
  for(i in seq_len(nsim)) {
    X <- datagen.rpoisppOnLines(lambda, S, ...)
    Y <- lpp(X, L)
    if(bugout) return(Y)
    result[[i]] <- Y
  }
  result <- as.solist(result)
  if(nsim > 0) names(result) <- paste("Simulation", 1:nsim)
  return(result)
}

runiflpp <- function(n, L, nsim=1, drop=TRUE) {
  verifyclass(L, "linnet")
  result <- vector(mode="list", length=nsim)
  S <- as.psp(L)
  bugout <- (nsim == 1) && drop
  for(i in seq_len(nsim)) {
    X <- datagen.runifpointOnLines(n, S)
    Y <- lpp(X, L)
    if(bugout) return(Y)
    result[[i]] <- Y
  }
  result <- as.solist(result)
  if(nsim > 0) names(result) <- paste("Simulation", 1:nsim)
  return(result)
}
