#
#      distanxD.R
#
#      $Revision: 1.5 $     $Date: 2013/04/25 06:37:43 $
#
#      Interpoint distances for multidimensional points
#
#      Methods for pairdist, nndist, nnwhich, crossdist
#

pairdist.ppx <- function(X, ...) {
  verifyclass(X, "ppx")
  # extract point coordinates
  coo <- as.matrix(coords(X, ...))
  n <- nrow(coo)
  if(n == 0)
    return(matrix(numeric(0), nrow=0, ncol=0))
  return(as.matrix(dist(coo)))
}

crossdist.ppx <- function(X, Y, ...) {
  verifyclass(X, "ppx")
  verifyclass(Y, "ppx")
  # extract point coordinates
  cooX <- as.matrix(coords(X, ...))
  cooY <- as.matrix(coords(Y, ...))
  nX <- nrow(cooX)
  nY <- nrow(cooY)
  if(ncol(cooX) != ncol(cooY))
    stop("X and Y have different dimensions (different numbers of coordinates)")
  if(nX == 0 || nY == 0)
    return(matrix(numeric(0), nrow=nX, ncol=nY))
  coo <- rbind(cooX, cooY)
  dis <- as.matrix(dist(coo))
  ans <- dis[1:nX, nX + (1:nY)]
  return(ans)
}

nndist.ppx <- function(X, ..., k=1) {
  verifyclass(X, "ppx")

  # extract point coordinates
  coo <- as.matrix(coords(X, ...))
  n <- nrow(coo)
  m <- ncol(coo)

  if(m == 0) {
    warning("nndist.ppx: Zero-dimensional coordinates: returning NA")
    if(length(k) == 1)
      return(rep.int(NA_real_, n))
    else
      return(matrix(NA_real_, n, length(k)))
  }
  
  # k can be a single integer or an integer vector
  if(length(k) == 0)
    stop("k is an empty vector")
  else if(length(k) == 1) {
    if(k != round(k) || k <= 0)
      stop("k is not a positive integer")
  } else {
    if(any(k != round(k)) || any(k <= 0))
      stop(paste("some entries of the vector",
           sQuote("k"), "are not positive integers"))
  }
  k <- as.integer(k)
  kmax <- max(k)

  # trivial cases
  if(n <= 1) {
    # empty pattern => return numeric(0)
    # or pattern with only 1 point => return Inf
    nnd <- matrix(Inf, nrow=n, ncol=kmax)
    nnd <- nnd[,k, drop=TRUE]
    return(nnd)
  }
  
  # number of neighbours that are well-defined
  kmaxcalc <- min(n-1, kmax)

  # calculate k-nn distances for k <= kmaxcalc
  
  if(kmaxcalc == 1) {
    # calculate nearest neighbour distance only
    nnd<-numeric(n)
    o <- fave.order(coo[,1])
    big <- sqrt(.Machine$double.xmax)
    DUP <- spatstat.options("dupC")
    Cout <- .C("nndMD",
               n= as.integer(n),
               m=as.integer(m),
               x= as.double(t(coo[o,])),
               nnd= as.double(nnd),
               as.double(big),
               DUP=DUP)
#               PACKAGE="spatstat")
    nnd[o] <- Cout$nnd
  } else {
    # case kmaxcalc > 1
    nnd<-numeric(n * kmaxcalc)
    o <- fave.order(coo[,1])
    big <- sqrt(.Machine$double.xmax)
    DUP <- spatstat.options("dupC")
    Cout <- .C("knndMD",
               n    = as.integer(n),
               m    = as.integer(m),
               kmax = as.integer(kmaxcalc),
               x    = as.double(t(coo[o,])),
               nnd  = as.double(nnd),
               huge = as.double(big),
               DUP=DUP)
#               PACKAGE="spatstat")
    nnd <- matrix(nnd, nrow=n, ncol=kmaxcalc)
    nnd[o, ] <- matrix(Cout$nnd, nrow=n, ncol=kmaxcalc, byrow=TRUE)
  }

  # post-processing
  if(kmax > kmaxcalc) {
    # add columns of Inf's
    infs <- matrix(as.numeric(Inf), nrow=n, ncol=kmax-kmaxcalc)
    nnd <- cbind(nnd, infs)
  }

  if(length(k) < kmax) {
    # select only the specified columns
    nnd <- nnd[, k, drop=TRUE]
  }
  
  return(nnd)
}

nnwhich.ppx <- function(X, ..., k=1) {
  verifyclass(X, "ppx")
  # k can be a single integer or an integer vector
  if(length(k) == 0)
    stop("k is an empty vector")
  else if(length(k) == 1) {
    if(k != round(k) || k <= 0)
      stop("k is not a positive integer")
  } else {
    if(any(k != round(k)) || any(k <= 0))
      stop(paste("some entries of the vector",
           sQuote("k"), "are not positive integers"))
  }
  k <- as.integer(k)
  kmax <- max(k)

  # extract point coordinates
  coo <- coords(X, ...)
  n <- nrow(coo)
  m <- ncol(coo)
  
  if(m == 0) {
    warning("nnwhich.ppx: Zero-dimensional coordinates: returning NA")
    if(length(k) == 1)
      return(rep.int(NA_real_, n))
    else
      return(matrix(NA_real_, n, length(k)))
  }
  
  # special cases
  if(n <= 1) {
    # empty pattern => return integer(0)
    # or pattern with only 1 point => return NA
    nnw <- matrix(NA_integer_, nrow=n, ncol=kmax)
    nnw <- nnw[,k, drop=TRUE]
    return(nnw)
  }

  # number of neighbours that are well-defined
  kmaxcalc <- min(n-1, kmax)

  # identify k-nn for k <= kmaxcalc

  if(kmaxcalc == 1) {
    # identify nearest neighbour only
    nnw <- integer(n)
    o <- fave.order(coo[,1])
    big <- sqrt(.Machine$double.xmax)
    DUP <- spatstat.options("dupC")
    Cout <- .C("nnwMD",
               n = as.integer(n),
               m = as.integer(m),
               x = as.double(t(coo[o,])),
               nnd = as.double(numeric(n)),
               nnwhich = as.integer(nnw),
               huge = as.double(big),
               DUP=DUP)
#               PACKAGE="spatstat")
    witch <- Cout$nnwhich
    if(any(witch <= 0))
      stop("Internal error: non-positive index returned from C code")
    if(any(witch > n))
      stop("Internal error: index returned from C code exceeds n")
    nnw[o] <- o[witch]
  } else {
    # case kmaxcalc > 1
    nnw <- matrix(integer(n * kmaxcalc), nrow=n, ncol=kmaxcalc)
    o <- fave.order(coo[,1])
    big <- sqrt(.Machine$double.xmax)
    DUP <- spatstat.options("dupC")
    Cout <- .C("knnwMD",
               n = as.integer(n),
               m = as.integer(m),
               kmax = as.integer(kmaxcalc),
               x = as.double(t(coo[o,])),
               nnd = as.double(numeric(n * kmaxcalc)),
               nnwhich = as.integer(nnw),
               huge = as.double(big),
               DUP=DUP)
#               PACKAGE="spatstat")
    witch <- Cout$nnwhich
    witch <- matrix(witch, nrow=n, ncol=kmaxcalc, byrow=TRUE)
    if(any(witch <= 0))
      stop("Internal error: non-positive index returned from C code")
    if(any(witch > n))
      stop("Internal error: index returned from C code exceeds n")
    # convert back to original ordering
    nnw[o,] <- matrix(o[witch], nrow=n, ncol=kmaxcalc)
  }
  
  # post-processing
  if(kmax > kmaxcalc) {
    # add columns of NA's
    nas <- matrix(NA_integer_, nrow=n, ncol=kmax-kmaxcalc)
    nnw <- cbind(nnw, nas)
  }

  if(length(k) < kmax) {
    # select only the specified columns
    nnw <- nnw[, k, drop=TRUE]
  }
  return(nnw)
}

