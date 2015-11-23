#
#  Perfect Simulation 
#
#  $Revision: 1.20 $ $Date: 2015/11/23 07:02:21 $
#
#  rStrauss
#  rHardcore
#  rStraussHard
#  rDiggleGratton
#  rDGS
#  rPenttinen

rStrauss <- function(beta, gamma=1, R=0, W=owin(), expand=TRUE,
                     nsim=1, drop=TRUE) {

  if(!missing(W)) 
    verifyclass(W, "owin")

  check.1.real(beta)
  check.1.real(gamma)
  check.1.real(R)

  check.finite(beta)
  check.finite(gamma)
  check.finite(R)
  
  stopifnot(beta > 0)
  stopifnot(gamma >= 0)
  stopifnot(gamma <= 1)
  stopifnot(R >= 0)

  runif(1)

  Wsim <- expandwinPerfect(W, expand, rmhexpand(distance=2*R))
  xrange <- Wsim$xrange
  yrange <- Wsim$yrange

  result <- vector(mode="list", length=nsim)

  for(i in 1:nsim) {
    storage.mode(beta) <- storage.mode(gamma) <- storage.mode(R) <- "double"
    storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
    z <- .Call("PerfectStrauss",
               beta,
               gamma,
               R,
               xrange,
               yrange)

    X <- z[[1]]
    Y <- z[[2]]
    nout <- z[[3]]
    times <- c(start=z[[4]], end=z[[5]])
    
    if(nout<0)
      stop("internal error: copying failed in PerfectStrauss")

    seqn <- seq_len(nout)
    P <- ppp(X[seqn], Y[seqn], window=Wsim, check=FALSE)
    if(attr(Wsim, "changed"))
      P <- P[W]
    attr(P, "times") <- times

    if(nsim == 1 && drop) return(P)
    
    result[[i]] <- P
  }
  result <- as.solist(result)
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}

#  Perfect Simulation of Hardcore process

rHardcore <- function(beta, R=0, W=owin(), expand=TRUE, nsim=1, drop=TRUE) {
  if(!missing(W)) 
    verifyclass(W, "owin")

  check.1.real(beta)
  check.1.real(R)

  check.finite(beta)
  check.finite(R)

  stopifnot(beta > 0)
  stopifnot(R    >= 0)

  runif(1)

  Wsim <- expandwinPerfect(W, expand, rmhexpand(distance=2*R))
  xrange <- Wsim$xrange
  yrange <- Wsim$yrange

  result <- vector(mode="list", length=nsim)

  for(i in 1:nsim) {
    storage.mode(beta) <- storage.mode(R) <- "double"
    storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
    z <- .Call("PerfectHardcore",
               beta,
               R,
               xrange,
               yrange)

    X <- z[[1]]
    Y <- z[[2]]
    nout <- z[[3]]
    
    if(nout<0)
      stop("internal error: copying failed in PerfectHardcore")

    seqn <- seq_len(nout)
    P <- ppp(X[seqn], Y[seqn], window=Wsim, check=FALSE)
    if(attr(Wsim, "changed"))
      P <- P[W]

    if(nsim == 1 && drop) return(P)
    result[[i]] <- P
  }
  result <- as.solist(result)
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}

#
#  Perfect simulation of hybrid Strauss-Hardcore
#        provided gamma <= 1
#

rStraussHard <- function(beta, gamma=1, R=0, H=0, W=owin(),
                         expand=TRUE, nsim=1, drop=TRUE) {
  if(!missing(W)) 
    verifyclass(W, "owin")

  check.1.real(beta)
  check.1.real(gamma)
  check.1.real(R)
  check.1.real(H)

  check.finite(beta)
  check.finite(gamma)
  check.finite(R)
  check.finite(H)
  
  stopifnot(beta > 0)
  stopifnot(gamma >= 0)
  if(gamma > 1)
    stop("Sorry, perfect simulation is only implemented for gamma <= 1")
  stopifnot(R >= 0)
  stopifnot(H >= 0)
  stopifnot(H <= R)

  runif(1)

  Wsim <- expandwinPerfect(W, expand, rmhexpand(distance=2*R))
  xrange <- Wsim$xrange
  yrange <- Wsim$yrange

  result <- vector(mode="list", length=nsim)

  for(i in 1:nsim) {
    storage.mode(beta) <- storage.mode(gamma) <-
      storage.mode(R) <- storage.mode(H) <- "double"
    storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
    z <- .Call("PerfectStraussHard",
               beta,
               gamma,
               R,
               H,
               xrange,
               yrange)

    X <- z[[1]]
    Y <- z[[2]]
    nout <- z[[3]]

    if(nout<0)
      stop("internal error: copying failed in PerfectStraussHard")

    seqn <- seq_len(nout)
    P <- ppp(X[seqn], Y[seqn], window=Wsim, check=FALSE)
    if(attr(Wsim, "changed"))
      P <- P[W]

    if(nsim == 1 && drop) return(P)
    result[[i]] <- P
  }
  result <- as.solist(result)
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}

#
#  Perfect Simulation of Diggle-Gratton process
#

rDiggleGratton <- function(beta, delta, rho, kappa=1, W=owin(),
                           expand=TRUE, nsim=1, drop=TRUE) {
  if(!missing(W)) 
    verifyclass(W, "owin")

  check.1.real(beta)
  check.1.real(delta)
  check.1.real(rho)
  check.1.real(kappa)

  check.finite(beta)
  check.finite(delta)
  check.finite(rho)
  check.finite(kappa)

  stopifnot(beta > 0)
  stopifnot(delta >= 0)
  stopifnot(rho   >= 0)
  stopifnot(delta <= rho)
  stopifnot(kappa >= 0)

  runif(1)

  Wsim <- expandwinPerfect(W, expand, rmhexpand(distance=2*rho))
  xrange <- Wsim$xrange
  yrange <- Wsim$yrange

  result <- vector(mode="list", length=nsim)

  for(i in 1:nsim) {
    storage.mode(beta) <- "double"
    storage.mode(delta) <- storage.mode(rho) <- storage.mode(kappa) <- "double"
    storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
    z <- .Call("PerfectDiggleGratton",
               beta,
               delta,
               rho,
               kappa,
               xrange,
               yrange)

    X <- z[[1]]
    Y <- z[[2]]
    nout <- z[[3]]

    if(nout<0)
      stop("internal error: copying failed in PerfectDiggleGratton")

    seqn <- seq_len(nout)
    P <- ppp(X[seqn], Y[seqn], window=Wsim, check=FALSE)
    if(attr(Wsim, "changed"))
      P <- P[W]

    if(nsim == 1 && drop) return(P)
    result[[i]] <- P
  }
  result <- as.solist(result)
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}


#
#  Perfect Simulation of Diggle-Gates-Stibbard process
#

rDGS <- function(beta, rho, W=owin(), expand=TRUE, nsim=1, drop=TRUE) {
  if(!missing(W)) 
    verifyclass(W, "owin")

  check.1.real(beta)
  check.1.real(rho)

  check.finite(beta)
  check.finite(rho)

  stopifnot(beta > 0)
  stopifnot(rho  >= 0)

  runif(1)

  Wsim <- expandwinPerfect(W, expand, rmhexpand(distance=2*rho))
  xrange <- Wsim$xrange
  yrange <- Wsim$yrange

  result <- vector(mode="list", length=nsim)

  for(i in 1:nsim) {
    storage.mode(beta) <- "double"
    storage.mode(rho) <- "double"
    storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
    z <- .Call("PerfectDGS",
               beta,
               rho,
               xrange,
               yrange)

    X <- z[[1]]
    Y <- z[[2]]
    nout <- z[[3]]
    
    if(nout<0)
      stop("internal error: copying failed in PerfectDGS")
    
    seqn <- seq_len(nout)
    P <- ppp(X[seqn], Y[seqn], window=Wsim, check=FALSE)
    if(attr(Wsim, "changed"))
      P <- P[W]

    if(nsim == 1 && drop) return(P)
    result[[i]] <- P
  }
  result <- as.solist(result)
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}


#
#  Perfect Simulation of Penttinen process
#

rPenttinen <- function(beta, gamma=1, R, W=owin(),
                       expand=TRUE, nsim=1, drop=TRUE) {
  if(!missing(W)) 
    verifyclass(W, "owin")

  check.1.real(beta)
  check.1.real(gamma)
  check.1.real(R)

  check.finite(beta)
  check.finite(gamma)
  check.finite(R)

  stopifnot(beta > 0)
  stopifnot(gamma >= 0)
  stopifnot(gamma <= 1)
  stopifnot(R >= 0)

  runif(1)

  Wsim <- expandwinPerfect(W, expand, rmhexpand(distance=2*R))
  xrange <- Wsim$xrange
  yrange <- Wsim$yrange

  result <- vector(mode="list", length=nsim)

  for(i in 1:nsim) {
    storage.mode(beta) <- storage.mode(gamma) <- storage.mode(R) <- "double"
    storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
    z <- .Call("PerfectPenttinen",
               beta,
               gamma,
               R,
               xrange,
               yrange)

    X <- z[[1]]
    Y <- z[[2]]
    nout <- z[[3]]
    
    if(nout<0)
      stop("internal error: copying failed in PerfectPenttinen")
    
    seqn <- seq_len(nout)
    P <- ppp(X[seqn], Y[seqn], window=Wsim, check=FALSE)
    if(attr(Wsim, "changed"))
      P <- P[W]

    if(nsim == 1 && drop) return(P)
    result[[i]] <- P
  }
  result <- as.solist(result)
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}


## .......  utilities .................................

expandwinPerfect <- function(W, expand, amount) {
  ## expand 'W' if expand=TRUE according to default 'amount'
  ## or expand 'W' using rmhexpand(expand)
  if(!is.logical(expand)) {
    amount <- rmhexpand(expand)
    expand <- TRUE
  }
  changed <- FALSE
  if(expand) {
    W <- expand.owin(W, amount)
    changed <- TRUE
  }
  if(!is.rectangle(W)) {
    W <- as.rectangle(W)
    changed <- TRUE
    warning(paste("Simulation will be performed in the containing rectangle",
                  "and clipped to the original window."),
            call.=FALSE)
  }
  attr(W, "changed") <- changed
  return(W)
}
