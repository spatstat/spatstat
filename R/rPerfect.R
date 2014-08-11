#
#  Perfect Simulation 
#
#  $Revision: 1.14 $ $Date: 2012/10/13 05:49:28 $
#
#  rStrauss
#  rHardcore
#  rStraussHard
#  rDiggleGratton
#  rDGS

rStrauss <- function(beta, gamma=1, R=0, W=owin()) {
  if(!missing(W)) {
    verifyclass(W, "owin")
    if(W$type != "rectangle")
      stop("W must be a rectangle")
  }

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

  nothing <- runif(1)

  xrange <- W$xrange
  yrange <- W$yrange
  storage.mode(beta) <- storage.mode(gamma) <- storage.mode(R) <- "double"
  storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
  z <- .Call("PerfectStrauss",
             beta,
             gamma,
             R,
             xrange,
             yrange,
             PACKAGE="spatstat")

  X <- z[[1]]
  Y <- z[[2]]
  nout <- z[[3]]
  times <- c(start=z[[4]], end=z[[5]])

  if(nout<0)
    stop("internal error: copying failed in PerfectStrauss")

  seqn <- seq_len(nout)
  P <- ppp(X[seqn], Y[seqn], window=W, check=FALSE)
  attr(P, "times") <- times
  return(P)
}

#  Perfect Simulation of Hardcore process

rHardcore <- function(beta, R=0, W=owin()) {
  if(!missing(W)) {
    verifyclass(W, "owin")
    if(W$type != "rectangle")
      stop("W must be a rectangle")
  }

  check.1.real(beta)
  check.1.real(R)

  check.finite(beta)
  check.finite(R)

  stopifnot(beta > 0)
  stopifnot(R    >= 0)

  nothing <- runif(1)

  xrange <- W$xrange
  yrange <- W$yrange
  storage.mode(beta) <- storage.mode(R) <- "double"
  storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
  z <- .Call("PerfectHardcore",
             beta,
             R,
             xrange,
             yrange,
             PACKAGE="spatstat")

  X <- z[[1]]
  Y <- z[[2]]
  nout <- z[[3]]

  if(nout<0)
    stop("internal error: copying failed in PerfectHardcore")

  seqn <- seq_len(nout)
  P <- ppp(X[seqn], Y[seqn], window=W, check=FALSE)
  return(P)
}

#
#  Perfect simulation of hybrid Strauss-Hardcore
#        provided gamma <= 1
#

rStraussHard <- function(beta, gamma=1, R=0, H=0, W=owin()) {
  if(!missing(W)) {
    verifyclass(W, "owin")
    if(W$type != "rectangle")
      stop("W must be a rectangle")
  }

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

  nothing <- runif(1)

  xrange <- W$xrange
  yrange <- W$yrange
  storage.mode(beta) <- storage.mode(gamma) <-
    storage.mode(R) <- storage.mode(H) <- "double"
  storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
  z <- .Call("PerfectStraussHard",
             beta,
             gamma,
             R,
             H,
             xrange,
             yrange,
             PACKAGE="spatstat")

  X <- z[[1]]
  Y <- z[[2]]
  nout <- z[[3]]

  if(nout<0)
    stop("internal error: copying failed in PerfectStraussHard")

  seqn <- seq_len(nout)
  P <- ppp(X[seqn], Y[seqn], window=W, check=FALSE)
  return(P)
}


#
#  Perfect Simulation of Diggle-Gratton process
#

rDiggleGratton <- function(beta, delta, rho, kappa=1, W=owin()) {
  if(!missing(W)) {
    verifyclass(W, "owin")
    if(W$type != "rectangle")
      stop("W must be a rectangle")
  }

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

  nothing <- runif(1)

  xrange <- W$xrange
  yrange <- W$yrange
  storage.mode(beta) <- "double"
  storage.mode(delta) <- storage.mode(rho) <- storage.mode(kappa) <- "double"
  storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
  z <- .Call("PerfectDiggleGratton",
             beta,
             delta,
             rho,
             kappa,
             xrange,
             yrange,
             PACKAGE="spatstat")

  X <- z[[1]]
  Y <- z[[2]]
  nout <- z[[3]]

  if(nout<0)
    stop("internal error: copying failed in PerfectDiggleGratton")

  seqn <- seq_len(nout)
  P <- ppp(X[seqn], Y[seqn], window=W, check=FALSE)
  return(P)
}


#
#  Perfect Simulation of Diggle-Gates-Stibbard process
#

rDGS <- function(beta, rho, W=owin()) {
  if(!missing(W)) {
    verifyclass(W, "owin")
    if(W$type != "rectangle")
      stop("W must be a rectangle")
  }

  check.1.real(beta)
  check.1.real(rho)

  check.finite(beta)
  check.finite(rho)

  stopifnot(beta > 0)
  stopifnot(rho  >= 0)

  nothing <- runif(1)

  xrange <- W$xrange
  yrange <- W$yrange
  storage.mode(beta) <- "double"
  storage.mode(rho) <- "double"
  storage.mode(xrange) <- storage.mode(yrange) <- "double"
  
  z <- .Call("PerfectDGS",
             beta,
             rho,
             xrange,
             yrange,
             PACKAGE="spatstat")

  X <- z[[1]]
  Y <- z[[2]]
  nout <- z[[3]]

  if(nout<0)
    stop("internal error: copying failed in PerfectDGS")

  seqn <- seq_len(nout)
  P <- ppp(X[seqn], Y[seqn], window=W, check=FALSE)
  return(P)
}


