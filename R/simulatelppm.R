##
## simulatelppm.R
##
##  Simulation of lppm objects
##
##  $Revision: 1.7 $  $Date: 2018/05/12 16:14:05 $
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
  result <- simulationresult(result, nsim, drop)
  result <- timed(result, starttime=starttime)
  return(result)
}

