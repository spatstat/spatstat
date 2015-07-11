##
## simulatelppm.R
##
##  Simulation of lppm objects
##
##  $Revision: 1.6 $  $Date: 2015/07/11 08:19:26 $
##

simulate.lppm <- function(object, nsim=1, ...,
                          new.coef=NULL,
                          progress=(nsim > 1),
                          drop=FALSE) {
  starttime <- proc.time()
  if(!is.poisson(object$fit))
    stop("Simulation of non-Poisson models is not yet implemented")
  lambda <- predict(object, ..., new.coef=new.coef)
  lmax <- if(is.im(lambda)) max(lambda) else unlist(lapply(lambda, max))
  L <- as.linnet(object)
  result <- vector(mode="list", length=nsim)
  pstate <- list()
  for(i in seq_len(nsim)) {
    if(progress) pstate <- progressreport(i, nsim, state=pstate)
    result[[i]] <- rpoislpp(lambda, L, lmax=lmax)
  }
  if(nsim == 1 && drop) {
    result <- result[[1]]
  } else {
    result <- as.solist(result)
    if(nsim > 0)
      names(result) <- paste("Simulation", 1:nsim)
  }
  result <- timed(result, starttime=starttime)
  return(result)
}

