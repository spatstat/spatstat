#
#      distanxD.R
#
#      $Revision: 1.14 $     $Date: 2019/12/04 09:08:47 $
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
    if(length(k) == 1L)
      return(rep.int(NA_real_, n))
    else
      return(matrix(NA_real_, n, length(k)))
  }
  
  # k can be a single integer or an integer vector
  if(length(k) == 0)
    stop("k is an empty vector")
  else if(length(k) == 1L) {
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
  if(n <= 1L) {
    # empty pattern => return numeric(0)
    # or pattern with only 1 point => return Inf
    nnd <- matrix(Inf, nrow=n, ncol=kmax)
    nnd <- nnd[,k, drop=TRUE]
    return(nnd)
  }
  
  # number of neighbours that are well-defined
  kmaxcalc <- min(n-1L, kmax)

  # calculate k-nn distances for k <= kmaxcalc
  
  if(kmaxcalc == 1L) {
    # calculate nearest neighbour distance only
    nnd<-numeric(n)
    o <- fave.order(coo[,1L])
    big <- sqrt(.Machine$double.xmax)
    Cout <- .C("nndMD",
               n= as.integer(n),
               m=as.integer(m),
               x= as.double(t(coo[o,])),
               nnd= as.double(nnd),
               as.double(big),
               PACKAGE = "spatstat")
    nnd[o] <- Cout$nnd
  } else {
    # case kmaxcalc > 1
    nnd<-numeric(n * kmaxcalc)
    o <- fave.order(coo[,1L])
    big <- sqrt(.Machine$double.xmax)
    Cout <- .C("knndMD",
               n    = as.integer(n),
               m    = as.integer(m),
               kmax = as.integer(kmaxcalc),
               x    = as.double(t(coo[o,])),
               nnd  = as.double(nnd),
               huge = as.double(big),
               PACKAGE = "spatstat")
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
  else if(length(k) == 1L) {
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
    if(length(k) == 1L)
      return(rep.int(NA_real_, n))
    else
      return(matrix(NA_real_, n, length(k)))
  }
  
  # special cases
  if(n <= 1L) {
    # empty pattern => return integer(0)
    # or pattern with only 1 point => return NA
    nnw <- matrix(NA_integer_, nrow=n, ncol=kmax)
    nnw <- nnw[,k, drop=TRUE]
    return(nnw)
  }

  # number of neighbours that are well-defined
  kmaxcalc <- min(n-1L, kmax)

  # identify k-nn for k <= kmaxcalc

  if(kmaxcalc == 1L) {
    # identify nearest neighbour only
    nnw <- integer(n)
    o <- fave.order(coo[,1L])
    big <- sqrt(.Machine$double.xmax)
    Cout <- .C("nnwMD",
               n = as.integer(n),
               m = as.integer(m),
               x = as.double(t(coo[o,])),
               nnd = as.double(numeric(n)),
               nnwhich = as.integer(nnw),
               huge = as.double(big),
               PACKAGE = "spatstat")
    witch <- Cout$nnwhich
    if(any(witch <= 0))
      stop("Internal error: non-positive index returned from C code")
    if(any(witch > n))
      stop("Internal error: index returned from C code exceeds n")
    nnw[o] <- o[witch]
  } else {
    # case kmaxcalc > 1
    nnw <- matrix(integer(n * kmaxcalc), nrow=n, ncol=kmaxcalc)
    o <- fave.order(coo[,1L])
    big <- sqrt(.Machine$double.xmax)
    Cout <- .C("knnwMD",
               n = as.integer(n),
               m = as.integer(m),
               kmax = as.integer(kmaxcalc),
               x = as.double(t(coo[o,])),
               nnd = as.double(numeric(n * kmaxcalc)),
               nnwhich = as.integer(nnw),
               huge = as.double(big),
               PACKAGE = "spatstat")
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

nncross.ppx <- function(X, Y, iX=NULL, iY=NULL,
                        what = c("dist", "which"),
                        ...,
                        k = 1) {
  verifyclass(X, "ppx")
  verifyclass(Y, "ppx")
  what <- match.arg(what, several.ok=TRUE)
  want.dist  <- "dist" %in% what
  want.which <- "which" %in% what
  ## k can be a single integer or an integer vector
  if(length(k) == 0)
    stop("k is an empty vector")
  else if(length(k) == 1L) {
    if(k != round(k) || k <= 0)
      stop("k is not a positive integer")
  } else {
    if(any(k != round(k)) || any(k <= 0))
      stop(paste("some entries of the vector",
           sQuote("k"), "are not positive integers"))
  }
  k <- as.integer(k)
  kmax <- max(k)
  nk <- length(k)
  ## extract point coordinates
  cooX <- as.matrix(coords(X, ...))
  nX <- nrow(cooX)
  m  <- ncol(cooX)
  cooY <- as.matrix(coords(Y, ...))
  nY <- nrow(cooY)
  mY <- ncol(cooY)
  ## check dimensions
  if(mY != m)
    stop(paste("Point patterns have different spatial dimensions:",
               m, "!=", mY),
         call.=FALSE)
  if(m == 0) {
    warning("nncross.ppx: Zero-dimensional coordinates: returning NA")
    if(nk == 1L) {
      NND <- if(want.dist) rep.int(NA_real_,     nX) else 0
      NNW <- if(want.which) rep.int(NA_integer_, nX) else 0
    } else {
      NND <- if(want.dist) matrix(NA_real_,     nX, nk) else 0
      NNW <- if(want.which) matrix(NA_integer_, nX, nk) else 0
    }
    return(packupNNdata(NND, NNW, what, k))
  }
  ## trivial cases
  if(nX == 0L || nY == 0L) {
    NND <- matrix(Inf,         nrow=nX, ncol=nk)
    NNW <- matrix(NA_integer_, nrow=nX, ncol=nk)
    return(packupNNdata(NND, NNW, what, k))
  }

  ## exclusion arguments
  if(is.null(iX) != is.null(iY))
    stop("If one of iX, iY is given, then both must be given")
  exclude <- (!is.null(iX) || !is.null(iY))
  if(exclude) {
    stopifnot(is.integer(iX) && is.integer(iY))
    if(length(iX) != nX)
      stop("length of iX does not match the number of points in X")
    if(length(iY) != nY)
      stop("length of iY does not match the number of points in Y")
  }

  ## number of neighbours that are well-defined
  kmaxcalc <- min(nY, kmax)

  ## find k-nearest neighbours for k <= kmaxcalc
  oX <- fave.order(cooX[,1L])
  oY <- fave.order(cooY[,1L])
  big <- sqrt(.Machine$double.xmax)

  if(kmaxcalc == 1L) {
    ## find nearest neighbour only
    nnd <- numeric(nX)
    nnw <- integer(nX)
    if(!exclude) {
      Cout <- .C("nnXwMD",
                 m  =as.integer(m),
                 n1 = as.integer(nX),
                 x1 = as.double(t(cooX[oX,])),
                 n2 = as.integer(nY),
                 x2 = as.double(t(cooY[oY,])),
                 nnd = as.double(nnd),
                 nnwhich = as.integer(nnw),
                 as.double(big),
                 PACKAGE = "spatstat")
    } else {
      Cout <- .C("nnXxMD",
                 m  =as.integer(m),
                 n1 = as.integer(nX),
                 x1 = as.double(t(cooX[oX,])),
                 id1 = as.integer(iX[oX]),
                 n2 = as.integer(nY),
                 x2 = as.double(t(cooY[oY,])),
                 id2 = as.integer(iY[oY]),
                 nnd = as.double(nnd),
                 nnwhich = as.integer(nnw),
                 as.double(big),
                 PACKAGE = "spatstat")
    }
    if(want.dist)
      nnd[oX] <- Cout$nnd
    if(want.which) {
      witch <- Cout$nnwhich
      if(any(witch <= 0))
        stop("Internal error: non-positive index returned from C code")
      if(any(witch > nY))
        stop("Internal error: index returned from C code exceeds npoints(Y)")
      nnw[oX] <- oY[witch]
    }
  } else {
    ## k-nearest
    nnd <- matrix(0,  nX, kmaxcalc)
    nnw <- matrix(0L, nX, kmaxcalc)
    if(!exclude) {
      Cout <- .C("knnXwMD",
                 m       = as.integer(m),
                 n1      = as.integer(nX),
                 x1      = as.double(t(cooX[oX,])),
                 n2      = as.integer(nY),
                 x2      = as.double(t(cooY[oY,])),
                 kmax    = as.integer(kmaxcalc),
                 nnd     = as.double(nnd),
                 nnwhich = as.integer(nnw),
                 huge    = as.double(big),
                 PACKAGE = "spatstat")
    } else {
      Cout <- .C("knnXxMD",
                 m       = as.integer(m),
                 n1      = as.integer(nX),
                 x1      = as.double(t(cooX[oX,])),
                 id1     = as.integer(iX[oX]),
                 n2      = as.integer(nY),
                 x2      = as.double(t(cooY[oY,])),
                 id2     = as.integer(iY[oY]),
                 kmax    = as.integer(kmaxcalc),
                 nnd     = as.double(nnd),
                 nnwhich = as.integer(nnw),
                 huge    = as.double(big),
                 PACKAGE = "spatstat")
    }
    dust <- Cout$nnd
    witch <- Cout$nnwhich
    if(any(notfound <- (witch <= 0 | witch > nY))) {
      dust[notfound] <- Inf
      witch[notfound] <- NA
    }
    nnd[oX, ] <- matrix(dust, nrow=nX, ncol=kmaxcalc, byrow=TRUE)
    nnw[oX, ] <- matrix(oY[witch], nrow=nX, ncol=kmaxcalc, byrow=TRUE)
  }

  ## post-processing
  if(kmax > kmaxcalc) {
    ## add columns of Inf's/NA's
    if(want.dist) {
      infs <- matrix(as.numeric(Inf), nrow=nX, ncol=kmax-kmaxcalc)
      nnd <- cbind(nnd, infs)
    }
    if(want.which) {
      nas <- matrix(NA_integer_,      nrow=nX, ncol=kmax-kmaxcalc)
      nnw <- cbind(nnw, nas)
    }
  }
  if(length(k) < kmax) {
    ## select only the specified columns
    if(want.dist) nnd <- nnd[, k, drop=TRUE]
    if(want.which) nnw <- nnw[, k, drop=TRUE]
  }
  return(packupNNdata(nnd, nnw, what, k))  
}


packupNNdata <- function(NND, NNW, what, k) {
  result <- as.data.frame(list(dist=NND, which=NNW)[what])
  if(max(k) > 1L) {
    colnames(result) <- c(if("dist" %in% what) paste0("dist.", k) else NULL,
                          if("which" %in% what) paste0("which.",k) else NULL)
  }
  if(ncol(result) == 1L)
    result <- result[, , drop=TRUE]
  return(result)
}
