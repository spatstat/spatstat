
#
#      distances.R
#
#      $Revision: 1.38 $     $Date: 2012/03/18 09:47:52 $
#
#
#      Interpoint distances
#
#

pairdist <- function(X, ...) {
  UseMethod("pairdist")
}

pairdist.ppp <- function(X, ..., periodic=FALSE, method="C", squared=FALSE) {
  verifyclass(X, "ppp")
  if(!periodic)
    return(pairdist.default(X$x, X$y, method=method, squared=squared))
  # periodic case
  W <- X$window
  if(W$type != "rectangle")
    stop(paste("periodic edge correction can't be applied",
               "in a non-rectangular window"))
  wide <- diff(W$xrange)
  high <- diff(W$yrange)
  return(pairdist.default(X$x, X$y, period=c(wide,high),
                          method=method, squared=squared))
}


pairdist.default <-
  function(X, Y=NULL, ..., period=NULL, method="C", squared=FALSE)
{
  xy <- xy.coords(X,Y)[c("x","y")]
  x <- xy$x
  y <- xy$y

  n <- length(x)
  if(length(y) != n)
    stop("lengths of x and y do not match")

  # special cases
  if(n == 0)
    return(matrix(numeric(0), nrow=0, ncol=0))
  else if(n == 1)
    return(matrix(0,nrow=1,ncol=1))

  if((periodic<- !is.null(period))) {
    stopifnot(is.numeric(period))
    stopifnot(length(period) == 2 || length(period) == 1)
    stopifnot(all(period > 0))
    if(length(period) == 1) period <- rep(period, 2)
    wide <- period[1]
    high <- period[2]
  }

  switch(method,
         interpreted={
           xx <- matrix(rep(x, n), nrow = n)
           yy <- matrix(rep(y, n), nrow = n)
           if(!periodic) {
             d2 <- (xx - t(xx))^2 + (yy - t(yy))^2
           } else {
             dx <- xx - t(xx)
             dy <- yy - t(yy)
             dx2 <- pmin(dx^2, (dx + wide)^2, (dx - wide)^2)
             dy2 <- pmin(dy^2, (dy + high)^2, (dy - high)^2)
             d2 <- dx2 + dy2
           }
           if(squared)
             dout <- d2
           else
             dout <- sqrt(d2)
         },
         C={
           d <- numeric( n * n)
           DUP <- spatstat.options("dupC")
           if(!periodic) {
             z<- .C(if(squared) "pair2dist" else "pairdist",
                    n = as.integer(n),
                    x= as.double(x), y= as.double(y), d= as.double(d),
                    DUP=DUP,
                    PACKAGE="spatstat")
           } else {
             z<- .C(if(squared) "pairP2dist" else "pairPdist",
                    n = as.integer(n),
                    x= as.double(x), y= as.double(y),
                    xwidth=as.double(wide), yheight=as.double(high),
                    d= as.double(d),
                    DUP=DUP,
                    PACKAGE="spatstat")
           }
           dout <- matrix(z$d, nrow=n, ncol=n)
         },
         stop(paste("Unrecognised method", sQuote(method)))
       )
  return(dout)
}

nndist <- function(X, ...) {
  UseMethod("nndist")
}

nndist.ppp <- function(X, ..., k=1, method="C") {
  verifyclass(X, "ppp")
  if((narg <- length(list(...))) > 0) 
    warning(paste(narg, "unrecognised",
                  ngettext(narg, "argument was", "arguments were"),
                  "ignored"))
  return(nndist.default(X$x, X$y, k=k, method=method))
}

nndist.default <-
  function(X, Y=NULL, ..., k=1, method="C")
{
	#  computes the vector of nearest-neighbour distances 
	#  for the pattern of points (x[i],y[i])
	#
  xy <- xy.coords(X,Y)[c("x","y")]
  x <- xy$x
  y <- xy$y

  # validate
  n <- length(x)
  if(length(y) != n)
    stop("lengths of x and y do not match")
  
  # other arguments ignored
  if((narg <- length(list(...))) > 0) 
    warning(paste(narg, "unrecognised",
                  ngettext(narg, "argument was", "arguments were"),
                  "ignored"))
  
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
    switch(method,
         interpreted={
           #  matrix of squared distances between all pairs of points
           sq <- function(a, b) { (a-b)^2 }
           squd <-  outer(x, x, sq) + outer(y, y, sq)
           #  reset diagonal to a large value so it is excluded from minimum
           diag(squd) <- Inf
           #  nearest neighbour distances
           nnd <- sqrt(apply(squd,1,min))
         },
         C={
           nnd<-numeric(n)
           o <- fave.order(y)
           big <- sqrt(.Machine$double.xmax)
           DUP <- spatstat.options("dupC")
           z<- .C("nndistsort",
                  n= as.integer(n),
                  x= as.double(x[o]), y= as.double(y[o]), nnd= as.double(nnd),
                  as.double(big),
                  DUP=DUP,
                  PACKAGE="spatstat")
           nnd[o] <- z$nnd
         },
         stop(paste("Unrecognised method", sQuote(method)))
         )
  } else {
    # case kmaxcalc > 1
    switch(method,
           interpreted={
             if(n <= 1000) {
               # form n x n matrix of squared distances
               D2 <- pairdist.default(x, y, method=method, squared=TRUE)
               # find k'th smallest squared distance
               diag(D2) <- Inf
               NND2 <- t(apply(D2, 1, sort))[, 1:kmaxcalc]
               nnd <- sqrt(NND2)
             } else {
               # avoid creating huge matrix
               # handle one row of D at a time
               NND2 <- matrix(numeric(n * kmaxcalc), nrow=n, ncol=kmaxcalc)
               for(i in seq_len(n)) {
                 D2i <- (x - x[i])^2 + (y - y[i])^2
                 D2i[i] <- Inf
                 NND2[i,] <- sort(D2i)[1:kmaxcalc]
               }
               nnd <- sqrt(NND2)
             }
           },
           C={
             nnd<-numeric(n * kmaxcalc)
             o <- fave.order(y)
             big <- sqrt(.Machine$double.xmax)
             DUP <- spatstat.options("dupC")
             z<- .C("knndsort",
                    n    = as.integer(n),
                    kmax = as.integer(kmaxcalc),
                    x    = as.double(x[o]),
                    y    = as.double(y[o]),
                    nnd  = as.double(nnd),
                    huge = as.double(big),
                    DUP=DUP,
                    PACKAGE="spatstat")
             nnd <- matrix(nnd, nrow=n, ncol=kmaxcalc)
             nnd[o, ] <- matrix(z$nnd, nrow=n, ncol=kmaxcalc, byrow=TRUE)
           },
           stop(paste("Unrecognised method", sQuote(method)))
           )
  }

  # post-processing
  if(kmax > kmaxcalc) {
    # add columns of Inf
    nas <- matrix(Inf, nrow=n, ncol=kmax-kmaxcalc)
    nnd <- cbind(nnd, nas)
  }

  if(length(k) < kmax) {
    # select only the specified columns
    nnd <- nnd[, k, drop=TRUE]
  }
  
  return(nnd)
}

nnwhich <- function(X, ...) {
  UseMethod("nnwhich")
}

nnwhich.ppp <- function(X, ..., k=1, method="C") {
  verifyclass(X, "ppp")
  if((narg <- length(list(...))) > 0) 
    warning(paste(narg, "unrecognised",
                  ngettext(narg, "argument was", "arguments were"),
                  "ignored"))
  return(nnwhich.default(X$x, X$y, k=k, method=method))
}

nnwhich.default <-
  function(X, Y=NULL, ..., k=1, method="C")
{
	#  identifies nearest neighbour of each point in
	#  the pattern of points (x[i],y[i])
	#
  xy <- xy.coords(X,Y)[c("x","y")]
  x <- xy$x
  y <- xy$y

  # validate
  n <- length(x)
  if(length(y) != n)
    stop("lengths of x and y do not match")
  
  # other arguments ignored
  if((narg <- length(list(...))) > 0) 
    warning(paste(narg, "unrecognised",
                  ngettext(narg, "argument was", "arguments were"),
                  "ignored"))
  
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

  # special cases
  if(n <= 1) {
    # empty pattern => return integer(0)
    # or pattern with only 1 point => return NA
    nnw <- matrix(as.integer(NA), nrow=n, ncol=kmax)
    nnw <- nnw[,k, drop=TRUE]
    return(nnw)
  }

  # number of neighbours that are well-defined
  kmaxcalc <- min(n-1, kmax)

  # identify k-nn for k <= kmaxcalc

  if(kmaxcalc == 1) {
    # identify nearest neighbour only
    switch(method,
           interpreted={
             #  matrix of squared distances between all pairs of points
             sq <- function(a, b) { (a-b)^2 }
             squd <-  outer(x, x, sq) + outer(y, y, sq)
             #  reset diagonal to a large value so it is excluded from minimum
             diag(squd) <- Inf
             #  nearest neighbours
             nnw <- apply(squd,1,which.min)
           },
           C={
             nnw <- integer(n)
             o <- fave.order(y)
             big <- sqrt(.Machine$double.xmax)
             DUP <- spatstat.options("dupC")
             z<- .C("nnwhichsort",
                    n = as.integer(n),
                    x = as.double(x[o]),
                    y = as.double(y[o]),
                    nnwhich = as.integer(nnw),
                    huge = as.double(big),
                    DUP=DUP,
                    PACKAGE="spatstat")
             witch <- z$nnwhich # sic 
             if(any(witch <= 0))
               stop("Internal error: non-positive index returned from C code")
             if(any(witch > n))
               stop("Internal error: index returned from C code exceeds n")
             nnw[o] <- o[witch]
           },
           stop(paste("Unrecognised method", sQuote(method)))
           )
  } else {
    # case kmaxcalc > 1
    switch(method,
           interpreted={
             if(n <= 1000) {
               # form n x n matrix of squared distances
               D2 <- pairdist.default(x, y, method=method, squared=TRUE)
               # find k'th smallest squared distance
               diag(D2) <- Inf
               nnw <- t(apply(D2, 1, fave.order))[, 1:kmaxcalc]
             } else {
               # avoid creating huge matrix
               # handle one row of D at a time
               nnw <- matrix(as.integer(NA), nrow=n, ncol=kmaxcalc)
               for(i in seq_len(n)) {
                 D2i <- (x - x[i])^2 + (y - y[i])^2
                 D2i[i] <- Inf
                 nnw[i,] <- fave.order(D2i)[1:kmaxcalc]
               }      
             }
           },
           C={
             nnw <- matrix(integer(n * kmaxcalc), nrow=n, ncol=kmaxcalc)
             o <- fave.order(y)
             big <- sqrt(.Machine$double.xmax)
             DUP <- spatstat.options("dupC")
             z<- .C("knnsort",
                    n = as.integer(n),
                    kmax = as.integer(kmaxcalc),
                    x = as.double(x[o]),
                    y = as.double(y[o]),
                    nnd = as.double(numeric(n * kmaxcalc)),
                    nnwhich = as.integer(nnw),
                    huge = as.double(big),
                    DUP=DUP,
                    PACKAGE="spatstat")
             witch <- z$nnwhich # sic
             witch <- matrix(witch, nrow=n, ncol=kmaxcalc, byrow=TRUE)
             if(any(witch <= 0))
               stop("Internal error: non-positive index returned from C code")
             if(any(witch > n))
               stop("Internal error: index returned from C code exceeds n")
             # convert back to original ordering
             nnw[o,] <- matrix(o[witch], nrow=n, ncol=kmaxcalc)
           },
           stop(paste("Unrecognised method", sQuote(method)))
           )
  }
  
  # post-processing
  if(kmax > kmaxcalc) {
    # add columns of NA's
    nas <- matrix(as.numeric(NA), nrow=n, ncol=kmax-kmaxcalc)
    nnw <- cbind(nnw, nas)
  }

  if(length(k) < kmax) {
    # select only the specified columns
    nnw <- nnw[, k, drop=TRUE]
  }
  return(nnw)
}

crossdist <- function(X, Y, ...) {
  UseMethod("crossdist")
}

crossdist.ppp <- function(X, Y, ..., periodic=FALSE, method="C") {
  verifyclass(X, "ppp")
  Y <- as.ppp(Y)
  if(!periodic)
    return(crossdist.default(X$x, X$y, Y$x, Y$y, method=method))
  # periodic case
  WX <- X$window
  WY <- Y$window
  if(WX$type != "rectangle" || WY$type != "rectangle")
    stop(paste("periodic edge correction can't be applied",
               "in non-rectangular windows"))
  if(!is.subset.owin(WX,WY) || !is.subset.owin(WY, WX))
    stop(paste("periodic edge correction is not implemented",
               "for the case when X and Y lie in different rectangles"))
  wide <- diff(WX$xrange)
  high <- diff(WX$yrange)
  return(crossdist.default(X$x, X$y, Y$x, Y$y,
                           period=c(wide,high), method=method))
}

crossdist.default <-
  function(X, Y, x2, y2, ..., period=NULL, method="C")
{
  x1 <- X
  y1 <- Y
  # returns matrix[i,j] = distance from (x1[i],y1[i]) to (x2[j],y2[j])
  if(length(x1) != length(y1))
    stop("lengths of x and y do not match")
  if(length(x2) != length(y2))
    stop("lengths of x2 and y2 do not match")
  n1 <- length(x1)
  n2 <- length(x2)
  if(n1 == 0 || n2 == 0)
    return(matrix(numeric(0), nrow=n1, ncol=n2))

  if((periodic<- !is.null(period))) {
    stopifnot(is.numeric(period))
    stopifnot(length(period) == 2 || length(period) == 1)
    stopifnot(all(period > 0))
    if(length(period) == 1) period <- rep(period, 2)
    wide <- period[1]
    high <- period[2]
  }

   switch(method,
         interpreted = {
                 X1 <- matrix(rep(x1, n2), ncol = n2)
                 Y1 <- matrix(rep(y1, n2), ncol = n2)
                 X2 <- matrix(rep(x2, n1), ncol = n1)
                 Y2 <- matrix(rep(y2, n1), ncol = n1)
                 if(!periodic) 
                   d <- sqrt((X1 - t(X2))^2 + (Y1 - t(Y2))^2)
                 else {
                   dx <- X1 - t(X2)
                   dy <- Y1 - t(Y2)
                   dx2 <- pmin(dx^2, (dx + wide)^2, (dx - wide)^2)
                   dy2 <- pmin(dy^2, (dy + high)^2, (dy - high)^2)
                   d <- sqrt(dx2 + dy2)
                 }
                 return(d)
               },
               C = {
                 DUP <- spatstat.options("dupC")
                 if(!periodic) {
                   z<- .C("crossdist",
                          nfrom = as.integer(n1),
                          xfrom = as.double(x1),
                          yfrom = as.double(y1),
                          nto = as.integer(n2),
                          xto = as.double(x2),
                          yto = as.double(y2),
                          d = as.double(matrix(0, nrow=n1, ncol=n2)),
                          DUP=DUP,
                          PACKAGE="spatstat")
                 } else {
                   z<- .C("crossPdist",
                          nfrom = as.integer(n1),
                          xfrom = as.double(x1),
                          yfrom = as.double(y1),
                          nto = as.integer(n2),
                          xto = as.double(x2),
                          yto = as.double(y2),
                          xwidth = as.double(wide),
                          yheight = as.double(high),
                          d = as.double(matrix(0, nrow=n1, ncol=n2)),
                          DUP=DUP,
                          PACKAGE="spatstat")
                 }
                 return(matrix(z$d, nrow=n1, ncol=n2))
               },
               stop(paste("Unrecognised method", method))
               )
      }

