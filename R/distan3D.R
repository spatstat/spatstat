#
#      distan3D.R
#
#      $Revision: 1.12 $     $Date: 2014/11/10 11:01:39 $
#
#      Interpoint distances for 3D points
#
#      Methods for pairdist, nndist, nnwhich, crossdist
#

pairdist.pp3 <- function(X, ..., periodic=FALSE, squared=FALSE) {
  verifyclass(X, "pp3")
  # extract point coordinates
  xyz <- coords(X)
  n <- nrow(xyz)
  x <- xyz$x
  y <- xyz$y
  z <- xyz$z
  #   
  # special cases
  if(n == 0)
    return(matrix(numeric(0), nrow=0, ncol=0))
  else if(n == 1L)
    return(matrix(0,nrow=1L,ncol=1L))
  #
  if(!periodic) {
    Cout <- .C("D3pairdist",
               n = as.integer(n),
               x = as.double(x),
               y = as.double(y),
               z = as.double(z),
               squared = as.integer(squared),
               d = as.double(numeric(n*n)))
  } else {
    b <- as.box3(X)
    wide <- diff(b$xrange)
    high <- diff(b$yrange)
    deep <- diff(b$zrange)
    Cout <- .C("D3pairPdist",
               n = as.integer(n),
               x = as.double(x),
               y = as.double(y),
               z = as.double(z),
               xwidth=as.double(wide),
               yheight=as.double(high),
               zdepth=as.double(deep),
               squared = as.integer(squared),
               d= as.double(numeric(n*n)))
  }
  dout <- matrix(Cout$d, nrow=n, ncol=n)
  return(dout)
}

nndist.pp3 <- function(X, ..., k=1) {
  verifyclass(X, "pp3")

  if((narg <- length(list(...))) > 0) 
    warning(paste(narg, "unrecognised",
                  ngettext(narg, "argument was", "arguments were"),
                  "ignored"))

  # extract point coordinates
  xyz <- coords(X)
  n <- nrow(xyz)
  x <- xyz$x
  y <- xyz$y
  z <- xyz$z
  
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
    o <- fave.order(z)
    big <- sqrt(.Machine$double.xmax)
    Cout <- .C("nnd3D",
               n= as.integer(n),
               x= as.double(x[o]),
               y= as.double(y[o]),
               z= as.double(z[o]),
               nnd= as.double(nnd),
               nnwhich = as.integer(integer(1L)),
               huge=as.double(big))
    nnd[o] <- Cout$nnd
  } else {
    # case kmaxcalc > 1
    nnd<-numeric(n * kmaxcalc)
    o <- fave.order(z)
    big <- sqrt(.Machine$double.xmax)
    Cout <- .C("knnd3D",
               n    = as.integer(n),
               kmax = as.integer(kmaxcalc),
               x    = as.double(x[o]),
               y    = as.double(y[o]),
               z    = as.double(z[o]),
               nnd  = as.double(nnd),
               nnwhich = as.integer(integer(1L)),
               huge = as.double(big))
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

nnwhich.pp3 <- function(X, ..., k=1) {
  verifyclass(X, "pp3")
  if((narg <- length(list(...))) > 0) 
    warning(paste(narg, "unrecognised",
                  ngettext(narg, "argument was", "arguments were"),
                  "ignored"))
  
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
  xyz <- coords(X)
  n <- nrow(xyz)
  x <- xyz$x
  y <- xyz$y
  z <- xyz$z
  
  # special cases
  if(n <= 1L) {
    # empty pattern => return integer(0)
    # or pattern with only 1 point => return NA
    nnw <- matrix(as.integer(NA), nrow=n, ncol=kmax)
    nnw <- nnw[,k, drop=TRUE]
    return(nnw)
  }

  # number of neighbours that are well-defined
  kmaxcalc <- min(n-1L, kmax)

  # identify k-nn for k <= kmaxcalc

  if(kmaxcalc == 1L) {
    # identify nearest neighbour only
    nnw <- integer(n)
    o <- fave.order(z)
    big <- sqrt(.Machine$double.xmax)
    Cout <- .C("nnw3D",
               n = as.integer(n),
               x = as.double(x[o]),
               y = as.double(y[o]),
               z = as.double(z[o]),
               nnd = as.double(numeric(1L)),
               nnwhich = as.integer(nnw),
               huge = as.double(big))
    # [sic] Conversion from C to R indexing is done in C code.
    witch <- Cout$nnwhich
    if(any(witch <= 0))
      stop("Internal error: illegal index returned from C code")
    if(any(witch > n))
      stop("Internal error: index returned from C code exceeds n")
    nnw[o] <- o[witch]
  } else {
    # case kmaxcalc > 1
    nnw <- matrix(integer(n * kmaxcalc), nrow=n, ncol=kmaxcalc)
    o <- fave.order(z)
    big <- sqrt(.Machine$double.xmax)
    Cout <- .C("knnw3D",
               n = as.integer(n),
               kmax = as.integer(kmaxcalc),
               x = as.double(x[o]),
               y = as.double(y[o]),
               z = as.double(z[o]),
               nnd = as.double(numeric(1L)),
               nnwhich = as.integer(nnw),
               huge = as.double(big))
    # [sic] Conversion from C to R indexing is done in C code.
    witch <- Cout$nnwhich 
    witch <- matrix(witch, nrow=n, ncol=kmaxcalc, byrow=TRUE)
    if(any(witch <= 0))
      stop("Internal error: illegal index returned from C code")
    if(any(witch > n))
      stop("Internal error: index returned from C code exceeds n")
    # convert back to original ordering
    nnw[o,] <- matrix(o[witch], nrow=n, ncol=kmaxcalc)
  }
  
  # post-processing
  if(kmax > kmaxcalc) {
    # add columns of NA's
    nas <- matrix(as.integer(NA), nrow=n, ncol=kmax-kmaxcalc)
    nnw <- cbind(nnw, nas)
  }

  if(length(k) < kmax) {
    # select only the specified columns
    nnw <- nnw[, k, drop=TRUE]
  }
  return(nnw)
}

crossdist.pp3 <- function(X, Y, ..., periodic=FALSE, squared=FALSE) {
  verifyclass(X, "pp3")
  verifyclass(Y, "pp3")

  cX <- coords(X)
  cY <- coords(Y)
  nX <- nrow(cX)
  nY <- nrow(cY)

  if(nX == 0 || nY == 0)
    return(matrix(numeric(0), nrow=nX, ncol=nY))

  if(!periodic) {
    Cout <- .C("D3crossdist",
               nfrom = as.integer(nX),
               xfrom = as.double(cX$x),
               yfrom = as.double(cX$y),
               zfrom = as.double(cX$z),
               nto = as.integer(nY),
               xto = as.double(cY$x),
               yto = as.double(cY$y),
               zto = as.double(cY$z),
               squared = as.integer(squared),
               d = as.double(matrix(0, nrow=nX, ncol=nY)))
  } else {
    b <- as.box3(X)
    wide <- diff(b$xrange)
    high <- diff(b$yrange)
    deep <- diff(b$zrange)
    Cout <- .C("D3crossPdist",
               nfrom = as.integer(nX),
               xfrom = as.double(cX$x),
               yfrom = as.double(cX$y),
               zfrom = as.double(cX$z),
               nto = as.integer(nY),
               xto = as.double(cY$x),
               yto = as.double(cY$y),
               zto = as.double(cY$z),
               xwidth = as.double(wide),
               yheight = as.double(high),
               zheight = as.double(deep),
               squared = as.integer(squared),
               d = as.double(matrix(0, nrow=nX, ncol=nY)))
  }
  return(matrix(Cout$d, nrow=nX, ncol=nY))
}
  
