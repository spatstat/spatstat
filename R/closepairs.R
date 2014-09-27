#
# closepairs.R
#
#   $Revision: 1.29 $   $Date: 2014/06/20 07:53:39 $
#
#  simply extract the r-close pairs from a dataset
# 
#  Less memory-hungry for large patterns
#

closepairs <- function(X, rmax, ...) {
  UseMethod("closepairs")
}
  
closepairs.ppp <- function(X, rmax, ordered=TRUE,
                           what=c("all", "indices"), ...) {
  verifyclass(X, "ppp")
  what <- match.arg(what)
  stopifnot(is.numeric(rmax) && length(rmax) == 1)
  stopifnot(is.finite(rmax))
  stopifnot(rmax >= 0)
  npts <- npoints(X)
  null.answer <- switch(what,
                        all = {
                          list(i=integer(0),
                               j=integer(0),
                               xi=numeric(0),
                               yi=numeric(0),
                               xj=numeric(0),
                               yj=numeric(0),
                               dx=numeric(0),
                               dy=numeric(0),
                               d=numeric(0))
                        },
                        indices = {
                          list(i=integer(0),
                               j=integer(0))
                        })
   if(npts == 0)
     return(null.answer)
  # sort points by increasing x coordinate
  oo <- fave.order(X$x)
  Xsort <- X[oo]
  # First make an OVERESTIMATE of the number of pairs
  nsize <- ceiling(4 * pi * (npts^2) * (rmax^2)/area(Window(X)))
  nsize <- max(1024, nsize)
  # Now extract pairs
  if(spatstat.options("closepairs.newcode")) {
    # ------------------- use new faster code ---------------------
    x <- Xsort$x
    y <- Xsort$y
    r <- rmax
    ng <- nsize
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(r) <- "double"
    storage.mode(ng) <- "integer"
    switch(what,
           all = {
             z <- .Call("Vclosepairs",
                        xx=x, yy=y, rr=r, nguess=ng)
             if(length(z) != 9)
               stop("Internal error: incorrect format returned from Vclosepairs")
             i  <- z[[1]]  # NB no increment required
             j  <- z[[2]]
             xi <- z[[3]]
             yi <- z[[4]]
             xj <- z[[5]]
             yj <- z[[6]]
             dx <- z[[7]]
             dy <- z[[8]]
             d  <- z[[9]]
           },
           indices = {
             z <- .Call("VcloseIJpairs",
                        xx=x, yy=y, rr=r, nguess=ng)
             if(length(z) != 2)
               stop("Internal error: incorrect format returned from VcloseIJpairs")
             i  <- z[[1]]  # NB no increment required
             j  <- z[[2]]
           })

  } else {
    # ------------------- use older code --------------------------
    DUP <- spatstat.options("dupC")
    z <-
      .C("Fclosepairs",
         nxy=as.integer(npts),
         x=as.double(Xsort$x),
         y=as.double(Xsort$y),
         r=as.double(rmax),
         noutmax=as.integer(nsize), 
         nout=as.integer(integer(1)),
         iout=as.integer(integer(nsize)),
         jout=as.integer(integer(nsize)), 
         xiout=as.double(numeric(nsize)),
         yiout=as.double(numeric(nsize)),
         xjout=as.double(numeric(nsize)),
         yjout=as.double(numeric(nsize)),
         dxout=as.double(numeric(nsize)),
         dyout=as.double(numeric(nsize)),
         dout=as.double(numeric(nsize)),
         status=as.integer(integer(1)),
         DUP=DUP)

    if(z$status != 0) {
      # Guess was insufficient
      # Obtain an OVERCOUNT of the number of pairs
      # (to work around gcc bug #323)
      rmaxplus <- 1.25 * rmax
      nsize <- .C("paircount",
                  nxy=as.integer(npts),
                  x=as.double(Xsort$x),
                  y=as.double(Xsort$y),
                  rmaxi=as.double(rmaxplus),
                  count=as.integer(integer(1)),
                  DUP=DUP)$count
      if(nsize <= 0)
        return(null.answer)
      # add a bit more for safety
      nsize <- ceiling(1.1 * nsize) + 2 * npts
      # now extract points
      z <-
        .C("Fclosepairs",
           nxy=as.integer(npts),
           x=as.double(Xsort$x),
           y=as.double(Xsort$y),
           r=as.double(rmax),
           noutmax=as.integer(nsize), 
           nout=as.integer(integer(1)),
           iout=as.integer(integer(nsize)),
           jout=as.integer(integer(nsize)), 
           xiout=as.double(numeric(nsize)),
           yiout=as.double(numeric(nsize)),
           xjout=as.double(numeric(nsize)),
           yjout=as.double(numeric(nsize)),
           dxout=as.double(numeric(nsize)),
           dyout=as.double(numeric(nsize)),
           dout=as.double(numeric(nsize)),
           status=as.integer(integer(1)),
           DUP=DUP)
      if(z$status != 0)
        stop(paste("Internal error: C routine complains that insufficient space was allocated:", nsize))
    }
  # trim vectors to the length indicated
    npairs <- z$nout
    if(npairs <= 0)
      return(null.answer)
    actual <- seq_len(npairs)
    i  <- z$iout[actual]  # sic
    j  <- z$jout[actual]
    if(what == "all") {
      xi <- z$xiout[actual]
      yi <- z$yiout[actual]
      xj <- z$xjout[actual]
      yj <- z$yjout[actual]
      dx <- z$dxout[actual]
      dy <- z$dyout[actual]
      d <-  z$dout[actual]
    }
    # ------------------- end code switch ------------------------
  }
  
  # convert i,j indices to original sequence
  i <- oo[i]
  j <- oo[j]
  # are (i, j) and (j, i) equivalent?
  if(!ordered) {
    ok <- (i < j)
    i  <-  i[ok]
    j  <-  j[ok]
    if(what == "all") {
      xi <- xi[ok]
      yi <- yi[ok]
      xj <- xj[ok]
      yj <- yj[ok]
      dx <- dx[ok]
      dy <- dy[ok]
      d  <-  d[ok]
    }
  }
  # done
  switch(what,
         all = {
           answer <- list(i=i,
                          j=j,
                          xi=xi, 
                          yi=yi,
                          xj=xj,
                          yj=yj,
                          dx=dx,
                          dy=dy,
                          d=d)
         },
         indices = {
           answer <- list(i = i, j = j)
         })
  return(answer)
}

#######################

crosspairs <- function(X, Y, rmax, ...) {
  UseMethod("crosspairs")
}

crosspairs.ppp <- function(X, Y, rmax, what=c("all", "indices"), ...) {
  verifyclass(X, "ppp")
  verifyclass(Y, "ppp")
  what <- match.arg(what)
  stopifnot(is.numeric(rmax) && length(rmax) == 1 && rmax >= 0)
  null.answer <- switch(what,
                        all = {
                          list(i=integer(0),
                               j=integer(0),
                               xi=numeric(0),
                               yi=numeric(0),
                               xj=numeric(0),
                               yj=numeric(0),
                               dx=numeric(0),
                               dy=numeric(0),
                               d=numeric(0))
                        },
                        indices = {
                          list(i=integer(0),
                               j=integer(0))
                        })
  nX <- npoints(X)
  nY <- npoints(Y)
  if(nX == 0 || nY == 0) return(null.answer)
  # order patterns by increasing x coordinate
  ooX <- fave.order(X$x)
  Xsort <- X[ooX]
  ooY <- fave.order(Y$x)
  Ysort <- Y[ooY]
  if(spatstat.options("crosspairs.newcode")) {
    # ------------------- use new faster code ---------------------
    # First (over)estimate the number of pairs
    nsize <- ceiling(2 * pi * (rmax^2) * nX * nY/area(Window(Y)))
    nsize <- max(1024, nsize)
    # .Call
    Xx <- Xsort$x
    Xy <- Xsort$y
    Yx <- Ysort$x
    Yy <- Ysort$y
    r <- rmax
    ng <- nsize
    storage.mode(Xx) <- storage.mode(Xy) <- "double"
    storage.mode(Yx) <- storage.mode(Yy) <- "double"
    storage.mode(r) <- "double"
    storage.mode(ng) <- "integer"
    switch(what,
           all = {
             z <- .Call("Vcrosspairs",
                        xx1=Xx, yy1=Xy,
                        xx2=Yx, yy2=Yy,
                        rr=r, nguess=ng)
             if(length(z) != 9)
               stop("Internal error: incorrect format returned from Vcrosspairs")
             i  <- z[[1]]  # NB no increment required
             j  <- z[[2]]
             xi <- z[[3]]
             yi <- z[[4]]
             xj <- z[[5]]
             yj <- z[[6]]
             dx <- z[[7]]
             dy <- z[[8]]
             d  <- z[[9]]
           },
           indices = {
             z <- .Call("VcrossIJpairs",
                        xx1=Xx, yy1=Xy,
                        xx2=Yx, yy2=Yy,
                        rr=r, nguess=ng)
             if(length(z) != 2)
               stop("Internal error: incorrect format returned from VcrossIJpairs")
             i  <- z[[1]]  # NB no increment required
             j  <- z[[2]]
           })
           
  } else {
    # Older code 
    # obtain upper estimate of number of pairs
    # (to work around gcc bug 323)
    DUP <- spatstat.options("dupC")
    rmaxplus <- 1.25 * rmax
    nsize <- .C("crosscount",
                nn1=as.integer(X$n),
                x1=as.double(Xsort$x),
                y1=as.double(Xsort$y),
                nn2=as.integer(Ysort$n),
                x2=as.double(Ysort$x),
                y2=as.double(Ysort$y),
                rmaxi=as.double(rmaxplus),
                count=as.integer(integer(1)),
                DUP=DUP)$count
    if(nsize <= 0)
      return(null.answer)

    # allow slightly more space to work around gcc bug #323
    nsize <- ceiling(1.1 * nsize) + X$n + Y$n
    
    # now extract pairs
    z <-
      .C("Fcrosspairs",
         nn1=as.integer(X$n),
         x1=as.double(Xsort$x),
         y1=as.double(Xsort$y),
         nn2=as.integer(Y$n),
         x2=as.double(Ysort$x),
         y2=as.double(Ysort$y),
         r=as.double(rmax),
         noutmax=as.integer(nsize), 
         nout=as.integer(integer(1)),
         iout=as.integer(integer(nsize)),
         jout=as.integer(integer(nsize)), 
         xiout=as.double(numeric(nsize)),
         yiout=as.double(numeric(nsize)),
         xjout=as.double(numeric(nsize)),
         yjout=as.double(numeric(nsize)),
         dxout=as.double(numeric(nsize)),
         dyout=as.double(numeric(nsize)),
         dout=as.double(numeric(nsize)),
         status=as.integer(integer(1)),
         DUP=DUP)
    if(z$status != 0)
      stop(paste("Internal error: C routine complains that insufficient space was allocated:", nsize))
    # trim vectors to the length indicated
    npairs <- z$nout
    if(npairs <= 0)
      return(null.answer)
    actual <- seq_len(npairs)
    i  <- z$iout[actual] # sic
    j  <- z$jout[actual] 
    xi <- z$xiout[actual]
    yi <- z$yiout[actual]
    xj <- z$xjout[actual]
    yj <- z$yjout[actual]
    dx <- z$dxout[actual]
    dy <- z$dyout[actual]
    d <-  z$dout[actual]
  }
  # convert i,j indices to original sequences
  i <- ooX[i]
  j <- ooY[j]
  # done
  switch(what,
         all = {
           answer <- list(i=i,
                          j=j,
                          xi=xi, 
                          yi=yi,
                          xj=xj,
                          yj=yj,
                          dx=dx,
                          dy=dy,
                          d=d)
         },
         indices = {
           answer <- list(i=i, j=j)
         })
  return(answer)
}

closethresh <- function(X, R, S, ordered=TRUE) {
  # list all R-close pairs
  # and indicate which of them are S-close (S < R)
  # so that results are consistent with closepairs(X,S)
  verifyclass(X, "ppp")
  stopifnot(is.numeric(R) && length(R) == 1 && R >= 0)
  stopifnot(is.numeric(S) && length(S) == 1 && S >= 0)
  stopifnot(S < R)
  npts <- npoints(X)
   if(npts == 0)
     return(list(i=integer(0), j=integer(0), t=logical(0)))
  # sort points by increasing x coordinate
  oo <- fave.order(X$x)
  Xsort <- X[oo]
  # First make an OVERESTIMATE of the number of pairs
  nsize <- ceiling(4 * pi * (npts^2) * (R^2)/area(Window(X)))
  nsize <- max(1024, nsize)
  # Now extract pairs
  x <- Xsort$x
  y <- Xsort$y
  r <- R
  s <- S
  ng <- nsize
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
  storage.mode(r) <- "double"
  storage.mode(s) <- "double"
  storage.mode(ng) <- "integer"
  z <- .Call("Vclosethresh",
             xx=x, yy=y, rr=r, ss=s, nguess=ng)
  if(length(z) != 3)
    stop("Internal error: incorrect format returned from Vclosethresh")
  i  <- z[[1]]  # NB no increment required
  j  <- z[[2]]
  th <- as.logical(z[[3]])
  
  # convert i,j indices to original sequence
  i <- oo[i]
  j <- oo[j]
  # are (i, j) and (j, i) equivalent?
  if(!ordered) {
    ok <- (i < j)
    i  <-  i[ok]
    j  <-  j[ok]
    th <- th[ok]
  }
  # done
  return(list(i=i, j=j, th=th))
}

crosspairquad <- function(Q, rmax, what=c("all", "indices")) {
  # find all close pairs X[i], U[j]
  stopifnot(inherits(Q, "quad"))
  what <- match.arg(what)
  X <- Q$data
  D <- Q$dummy
  clX <- closepairs(X=X, rmax=rmax, what=what)
  clXD <- crosspairs(X=X, Y=D, rmax=rmax, what=what)
  # convert all indices to serial numbers in union.quad(Q)
  # assumes data are listed first
  clXD$j <- npoints(X) + clXD$j
  result <- list(rbind(as.data.frame(clX), as.data.frame(clXD)))
  return(result)
}

  
