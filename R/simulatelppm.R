##
## simulatelppm.R
##
##  Simulation of lppm objects
##
##  $Revision: 1.5 $  $Date: 2014/06/06 03:10:41 $
##

simulate.lppm <- function(object, nsim=1, ...,
                          new.coef=NULL,
                          progress=(nsim > 1)) {
  starttime <- proc.time()
  if(!is.poisson(object$fit))
    stop("Simulation of non-Poisson models is not yet implemented")
  lambda <- predict(object, ..., new.coef=new.coef)
  lmax <- if(is.im(lambda)) max(lambda) else unlist(lapply(lambda, max))
  L <- as.linnet(object)
  result <- vector(mode="list", length=nsim)
  for(i in seq_len(nsim)) {
    if(progress) progressreport(i, nsim)
    result[[i]] <- rpoislpp(lambda, L, lmax=lmax)
  }
  result <- as.listof(result)
  result <- timed(result, starttime=starttime)
  return(result)
}

