#'
#'   randompp3.R
#'
#'   $Revision: 1.1 $ $Date: 2020/11/30 11:43:50 $
#'

runifpoint3 <- function(n, domain=box3(), nsim=1, drop=TRUE) {
  domain <- as.box3(domain)
  result <- vector(mode="list", length=nsim)
  dd <- as.list(domain)[c("xrange", "yrange", "zrange")]
  for(i in 1:nsim) {
    x <- with(dd, runif(n, min=xrange[1], max=xrange[2]))
    y <- with(dd, runif(n, min=yrange[1], max=yrange[2]))
    z <- with(dd, runif(n, min=zrange[1], max=zrange[2]))
    result[[i]] <- pp3(x,y,z,domain)
  }
  if(drop && nsim == 1) return(result[[1]])
  result <- as.anylist(result)
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}

rpoispp3 <- function(lambda, domain=box3(), nsim=1, drop=TRUE) {
  domain <- as.box3(domain)
  v <- volume(domain)
  if(!(is.numeric(lambda) && length(lambda) == 1))
    stop("lambda must be a single numeric value")
  np <- rpois(nsim, lambda * v)
  dd <- as.list(domain)[c("xrange", "yrange", "zrange")]
  result <- vector(mode="list", length=nsim)
  for(i in 1:nsim) {
    ni <- np[i]
    x <- with(dd, runif(ni, min=xrange[1], max=xrange[2]))
    y <- with(dd, runif(ni, min=yrange[1], max=yrange[2]))
    z <- with(dd, runif(ni, min=zrange[1], max=zrange[2]))
    result[[i]] <- pp3(x,y,z,domain)
  }
  if(drop && nsim == 1) return(result[[1]])
  result <- as.anylist(result)
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}

