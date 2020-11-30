#'
#'   randomppx.R
#'
#'   $Revision: 1.1 $ $Date: 2020/11/30 11:44:46 $
#'

runifpointx <- function(n, domain, nsim=1, drop=TRUE) {
  check.1.integer(n)
  check.1.integer(nsim)
  stopifnot(inherits(domain, "boxx"))
  ra <- domain$ranges
  d <- length(ra)
  result <- vector(mode="list", length=nsim)
  for(i in 1:nsim) {
    if(n == 0) {
      coo <- matrix(numeric(0), nrow=0, ncol=d)
    } else {
      coo <- mapply(runif,
                    n=rep(n, d),
                    min=ra[1,],
                    max=ra[2,])
      if(!is.matrix(coo)) coo <- matrix(coo, ncol=d)
    }
    colnames(coo) <- colnames(ra)
    df <- as.data.frame(coo)
    result[[i]] <- ppx(df, domain, coord.type=rep("s", d))
  }
  if(nsim == 1 && drop)
    return(result[[1]])
  result <- as.anylist(result)
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}

rpoisppx <- function(lambda, domain, nsim=1, drop=TRUE) {
  stopifnot(inherits(domain, "boxx"))
  stopifnot(is.numeric(lambda) && length(lambda) == 1 && lambda >= 0)
  n <- rpois(nsim, lambda * volume.boxx(domain))
  result <- vector(mode="list", length=nsim)
  for(i in 1:nsim) 
    result[[i]] <- runifpointx(n[i], domain)
  if(nsim == 1 && drop)
    return(result[[1]])
  result <- as.anylist(result)
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}

