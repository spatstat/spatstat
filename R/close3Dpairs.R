#
# close3Dpairs.R
#
#   $Revision: 1.5 $   $Date: 2014/06/20 09:11:03 $
#
#  extract the r-close pairs from a 3D dataset
# 
#
closepairs.pp3 <- function(X, rmax, ordered=TRUE,
                           what=c("all", "indices"), ...) {
  verifyclass(X, "pp3")
  what <- match.arg(what)
  stopifnot(is.numeric(rmax) && length(rmax) == 1)
  stopifnot(is.finite(rmax))
  stopifnot(rmax >= 0)
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
    null.answer <- lapply(nama, function(e) numeric(0))
    return(null.answer)
  }
  # sort points by increasing x coordinate
  oo <- fave.order(coords(X)$x)
  Xsort <- X[oo]
  # First make an OVERESTIMATE of the number of pairs
  nsize <- ceiling(5 * pi * (npts^2) * (rmax^3)/volume(as.box3(X)))
  nsize <- max(1024, nsize)
  # Now extract pairs
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
                      xx=x, yy=y, zz=z, rr=r, nguess=ng)
              },
              indices = {
                .Call("close3IJpairs",
                      xx=x, yy=y, zz=z, rr=r, nguess=ng)
              })
  names(a) <- nama
  # convert i,j indices to original sequence
  a$i <- oo[a$i]
  a$j <- oo[a$j]
  # are (i, j) and (j, i) equivalent?
  if(!ordered) {
    a <- as.data.frame(a)
    a <- a[with(a, i < j), , drop=FALSE]
    a <- as.list(a)
  }
  return(a)
}

#######################

crosspairs.pp3 <- function(X, Y, rmax, what=c("all", "indices"), ...) {
  verifyclass(X, "pp3")
  verifyclass(Y, "pp3")
  what <- match.arg(what)
  stopifnot(is.numeric(rmax) && length(rmax) == 1 && rmax >= 0)
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
    null.answer <- lapply(nama, function(e) numeric(0))
    return(null.answer)
  }
  # order patterns by increasing x coordinate
  ooX <- fave.order(coords(X)$x)
  Xsort <- X[ooX]
  ooY <- fave.order(coords(Y)$x)
  Ysort <- Y[ooY]
  ## First (over)estimate the number of pairs
  nsize <- ceiling(3 * pi * (rmax^3) * nX * nY/volume(as.box3(Y)))
  nsize <- max(1024, nsize)
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
                      rr=r, nguess=ng)
              },
              indices = {
                .Call("cross3IJpairs",
                      xx1=Xx, yy1=Xy, zz1=Xz,
                      xx2=Yx, yy2=Yy, zz2=Yz,
                      rr=r, nguess=ng)
           })
  names(a) <- nama
  # convert i,j indices to original sequence
  a$i <- ooX[a$i]
  a$j <- ooY[a$j]
  return(a)
}

