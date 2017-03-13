#
# close3Dpairs.R
#
#   $Revision: 1.8 $   $Date: 2016/03/06 11:18:54 $
#
#  extract the r-close pairs from a 3D dataset
# 
#
closepairs.pp3 <- local({

  closepairs.pp3 <- function(X, rmax, twice=TRUE,
                             what=c("all", "indices"),
                             distinct=TRUE, neat=TRUE, ...) {
    verifyclass(X, "pp3")
    what <- match.arg(what)
    stopifnot(is.numeric(rmax) && length(rmax) == 1L)
    stopifnot(is.finite(rmax))
    stopifnot(rmax >= 0)
    ordered <- list(...)$ordered
    if(missing(twice) && !is.null(ordered)) {
      warning("Obsolete argument 'ordered' has been replaced by 'twice'")
      twice <- ordered
    }
    npts <- npoints(X)
    nama <- switch(what,
                   all = c("i", "j",
                           "xi", "yi", "zi",
                           "xj", "yj", "zj",
                           "dx", "dy", "dz",
                           "d"),
                   indices = c("i", "j"))
    names(nama) <- nama
    if(npts == 0) {
      null.answer <- lapply(nama, nuttink)
      return(null.answer)
    }
    ## sort points by increasing x coordinate
    oo <- fave.order(coords(X)$x)
    Xsort <- X[oo]
    ## First make an OVERESTIMATE of the number of pairs
    nsize <- ceiling(5 * pi * (npts^2) * (rmax^3)/volume(as.box3(X)))
    nsize <- max(1024, nsize)
    if(nsize > .Machine$integer.max) {
      warning(
        "Estimated number of close pairs exceeds maximum possible integer",
        call.=FALSE)
      nsize <- .Machine$integer.max
    }
    ## Now extract pairs
    XsortC <- coords(Xsort)
    x <- XsortC$x
    y <- XsortC$y
    z <- XsortC$z
    r <- rmax
    ng <- nsize
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(z) <- "double"
    storage.mode(r) <- "double"
    storage.mode(ng) <- "integer"
    ## go
    a <- switch(what,
                all = {
                  .Call("close3pairs",
                        xx=x, yy=y, zz=z, rr=r, nguess=ng,
                        PACKAGE = "spatstat")
                },
                indices = {
                  .Call("close3IJpairs",
                        xx=x, yy=y, zz=z, rr=r, nguess=ng,
                        PACKAGE = "spatstat")
                })
    names(a) <- nama
    ## convert i,j indices to original sequence
    a$i <- oo[a$i]
    a$j <- oo[a$j]
    ## handle options
    if(twice) {
      ## both (i, j) and (j, i) should be returned
      a <- as.data.frame(a)
      a <- as.list(rbind(a, swapdata(a, what)))
    } else if(neat) {
      ## enforce i < j
      swap <- with(a, (j < i))
      if(any(swap)) {
        a <- as.data.frame(a)
        a[swap,] <- swapdata(a[swap, ,drop=FALSE], what)
        a <- as.list(a)
      }
    }
    ## add pairs of identical points?
    if(!distinct) {
      ii <- seq_len(npts)
      xtra <- data.frame(i = ii, j=ii)
      if(what == "all") {
        coo <- coords(X)[, c("x","y","z")]
        zeroes <- rep(0, npts)
        xtra <- cbind(xtra, coo, coo, zeroes, zeroes, zeroes, zeroes)
      }
      a <- as.list(rbind(as.data.frame(a), xtra))
    }
    ## done
    return(a)
  }

  swapdata <- function(a, what) {
    switch(what,
           all = {
             with(a, data.frame(i  =  j,
                                j  =  i,
                                xi =  xj,
                                yi =  yj,
                                zi =  zj,
                                xj =  xi,
                                yj =  yi,
                                zj =  zi,
                                dx = -dx,
                                dy = -dy,
                                dz = -dz,
                                d  =  d))
           },
           indices = {
             with(a, data.frame(i=j,
                                j=i))
           })
  }
  
  nuttink <- function(x) numeric(0)

  closepairs.pp3
})

#######################

crosspairs.pp3 <- local({

  crosspairs.pp3 <- function(X, Y, rmax, what=c("all", "indices"), ...) {
    verifyclass(X, "pp3")
    verifyclass(Y, "pp3")
    what <- match.arg(what)
    stopifnot(is.numeric(rmax) && length(rmax) == 1L && rmax >= 0)
    nama <- switch(what,
                   all = c("i", "j",
                           "xi", "yi", "zi",
                           "xj", "yj", "zj",
                           "dx", "dy", "dz",
                           "d"),
                   indices = c("i", "j"))
    names(nama) <- nama
    nX <- npoints(X)
    nY <- npoints(Y)
    if(nX == 0 || nY == 0) {
      null.answer <- lapply(nama, nuttink)
      return(null.answer)
    }
    ## order patterns by increasing x coordinate
    ooX <- fave.order(coords(X)$x)
    Xsort <- X[ooX]
    ooY <- fave.order(coords(Y)$x)
    Ysort <- Y[ooY]
    ## First (over)estimate the number of pairs
    nsize <- ceiling(3 * pi * (rmax^3) * nX * nY/volume(as.box3(Y)))
    nsize <- max(1024, nsize)
    if(nsize > .Machine$integer.max) {
      warning(
        "Estimated number of close pairs exceeds maximum possible integer",
        call.=FALSE)
      nsize <- .Machine$integer.max
    }
    ## .Call
    XsortC <- coords(Xsort)
    YsortC <- coords(Ysort)
    Xx <- XsortC$x
    Xy <- XsortC$y
    Xz <- XsortC$z
    Yx <- YsortC$x
    Yy <- YsortC$y
    Yz <- YsortC$z
    r <- rmax
    ng <- nsize
    storage.mode(Xx) <- storage.mode(Xy) <- storage.mode(Xz) <- "double"
    storage.mode(Yx) <- storage.mode(Yy) <- storage.mode(Yz) <- "double"
    storage.mode(r) <- "double"
    storage.mode(ng) <- "integer"
    ## go
    a <- switch(what,
                all = {
                  .Call("cross3pairs",
                        xx1=Xx, yy1=Xy, zz1=Xz,
                        xx2=Yx, yy2=Yy, zz2=Yz,
                        rr=r, nguess=ng,
                        PACKAGE = "spatstat")
                },
                indices = {
                  .Call("cross3IJpairs",
                        xx1=Xx, yy1=Xy, zz1=Xz,
                        xx2=Yx, yy2=Yy, zz2=Yz,
                        rr=r, nguess=ng,
                        PACKAGE = "spatstat")
                })
    names(a) <- nama
    ## convert i,j indices to original sequence
    a$i <- ooX[a$i]
    a$j <- ooY[a$j]
    return(a)
  }

  nuttink <- function(x) numeric(0)
  
  crosspairs.pp3
})


