#
#      distances.R
#
#      $Revision: 1.44 $     $Date: 2013/11/03 05:20:13 $
#
#
#      Interpoint distances between pairs 
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
    if(length(period) == 1) period <- rep.int(period, 2)
    wide <- period[1]
    high <- period[2]
  }

  switch(method,
         interpreted={
           xx <- matrix(rep.int(x, n), nrow = n)
           yy <- matrix(rep.int(y, n), nrow = n)
           if(!periodic) {
             d2 <- (xx - t(xx))^2 + (yy - t(yy))^2
           } else {
             dx <- xx - t(xx)
             dy <- yy - t(yy)
             dx2 <- pmin.int(dx^2, (dx + wide)^2, (dx - wide)^2)
             dy2 <- pmin.int(dy^2, (dy + high)^2, (dy - high)^2)
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
               z<- .C("Cpairdist",
                      n = as.integer(n),
                      x= as.double(x),
                      y= as.double(y),
                      squared=as.integer(squared),
                      d= as.double(d),
                      DUP=DUP)
           } else {
             z <- .C("CpairPdist",
                     n = as.integer(n),
                     x= as.double(x),
                     y= as.double(y),
                     xwidth=as.double(wide),
                     yheight=as.double(high),
                     squared = as.integer(squared),
                     d= as.double(d),
                     DUP=DUP)
           }
           dout <- matrix(z$d, nrow=n, ncol=n)
         },
         stop(paste("Unrecognised method", sQuote(method)))
       )
  return(dout)
}


crossdist <- function(X, Y, ...) {
  UseMethod("crossdist")
}

crossdist.ppp <- function(X, Y, ...,
                          periodic=FALSE, method="C", squared=FALSE) {
  verifyclass(X, "ppp")
  Y <- as.ppp(Y)
  if(!periodic)
    return(crossdist.default(X$x, X$y, Y$x, Y$y,
                             method=method, squared=squared))
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
                           period=c(wide,high),
                           method=method, squared=squared))
}

crossdist.default <-
  function(X, Y, x2, y2, ..., period=NULL, method="C", squared=FALSE)
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
    if(length(period) == 1) period <- rep.int(period, 2)
    wide <- period[1]
    high <- period[2]
  }

   switch(method,
         interpreted = {
                 X1 <- matrix(rep.int(x1, n2), ncol = n2)
                 Y1 <- matrix(rep.int(y1, n2), ncol = n2)
                 X2 <- matrix(rep.int(x2, n1), ncol = n1)
                 Y2 <- matrix(rep.int(y2, n1), ncol = n1)
                 if(!periodic) 
                   d2 <- (X1 - t(X2))^2 + (Y1 - t(Y2))^2
                 else {
                   dx <- X1 - t(X2)
                   dy <- Y1 - t(Y2)
                   dx2 <- pmin.int(dx^2, (dx + wide)^2, (dx - wide)^2)
                   dy2 <- pmin.int(dy^2, (dy + high)^2, (dy - high)^2)
                   d2 <- dx2 + dy2
                 }
                 return(if(squared) d2 else sqrt(d2))
               },
               C = {
                 DUP <- spatstat.options("dupC")
                 if(!periodic) {
                   z<- .C("Ccrossdist",
                          nfrom = as.integer(n1),
                          xfrom = as.double(x1),
                          yfrom = as.double(y1),
                          nto = as.integer(n2),
                          xto = as.double(x2),
                          yto = as.double(y2),
                          squared = as.integer(squared),
                          d = as.double(matrix(0, nrow=n1, ncol=n2)),
                          DUP=DUP)
                 } else {
                   z<- .C("CcrossPdist",
                          nfrom = as.integer(n1),
                          xfrom = as.double(x1),
                          yfrom = as.double(y1),
                          nto = as.integer(n2),
                          xto = as.double(x2),
                          yto = as.double(y2),
                          xwidth = as.double(wide),
                          yheight = as.double(high),
                          squared = as.integer(squared),
                          d = as.double(matrix(0, nrow=n1, ncol=n2)),
                          DUP=DUP)
                 }
                 return(matrix(z$d, nrow=n1, ncol=n2))
               },
               stop(paste("Unrecognised method", method))
               )
      }

