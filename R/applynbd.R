# 	applynbd.R
#
#     $Revision: 1.15 $     $Date: 2013/04/25 06:37:43 $
#
#  applynbd()
# For each point, identify either
#	 - all points within distance R
#        - the closest N points  
#        - those points satisfying some constraint
# and apply the function FUN to them
#
#  markstat()
#      simple application of applynbd
#################################################################


applynbd <- function(X, FUN, N=NULL, R=NULL, criterion=NULL, exclude=FALSE, ...) {

  if(is.null(N) && is.null(R) && is.null(criterion)) 
    stop(paste("must specify at least one of the arguments",
               commasep(sQuote(c("N","R","criterion")))))
     
  X <- as.ppp(X)
  npts <- npoints(X)

  # compute matrix of pairwise distances
  dist <- pairdist(X)

  # compute row ranks (avoid ties)
  rankit <- function(x) {  u <- numeric(length(x)); u[fave.order(x)] <- seq_along(x); return(u) }
  drank <- t(apply(dist, 1, rankit)) - 1

  included <- matrix(TRUE, npts, npts)
  if(!is.null(R)) {
    # select points closer than R
    included <- included & (dist <= R)
  }
  if(!is.null(N)) {
    # select N closest points
    if(N < 1)
      stop("Value of N must be at least 1")
    if(exclude)
      included <- included & (drank <= N) 
    else
      included <- included & (drank <= N-1)
  }
  if(!is.null(criterion)) {
    # some funny criterion
    for(i in 1:npts) 
      included[i,] <- included[i,] & criterion(dist[i,], drank[i,])
  }
     
  if(exclude) 
    diag(included) <- FALSE

  # bind into an array
  a <- array(c(included, dist, drank, row(included)), dim=c(npts,npts,4))

  # what to do with a[i, , ]
  if(!is.marked(X)) 
    go <- function(ai, Z, fun, ...) { 
      which <- as.logical(ai[,1])
      distances <- ai[,2]
      dranks <- ai[,3]
      here <- ai[1,4]
      fun(Y=Z[which],
          current=c(x=Z$x[here], y=Z$y[here]),
          dists=distances[which], dranks=dranks[which],
          ...) 
    }
  else
    go <- function(ai, Z, fun, ...) { 
      which <- as.logical(ai[,1])
      distances <- ai[,2]
      dranks <- ai[,3]
      here <- ai[1,4]
      fun(Y=Z[which],
          current=Z[here],
          dists=distances[which], dranks=dranks[which],
          ...) 
    }
  
  # do it
  result <- apply(a, 1, go, Z=X, fun=FUN, ...)
  
  return(result)
}

markstat <- function(X, fun, N=NULL, R=NULL, ...) {
  verifyclass(X, "ppp")
  stopifnot(is.function(fun))
  statfun <- function(Y, current, dists, dranks, func, ...)
    { func(marks(Y, dfok=TRUE), ...) }
  applynbd(X, statfun, R=R, N=N, func=fun, ...)
}
