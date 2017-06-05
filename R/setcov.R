#
#
#     setcov.R
#
#     $Revision: 1.15 $ $Date: 2017/06/05 10:31:58 $
#
#    Compute the set covariance function of a window
#    or the (noncentred) spatial covariance function of an image
#

setcov <- function(W, V=W, ...) {
  W <- as.owin(W)
  # pixel approximation
  mW <- as.mask(W, ...)
  Z <- as.im(mW, na.replace=0)
  if(missing(V)) 
    return(imcov(Z))
  # cross-covariance
  V <- as.owin(V)
  mV <- as.mask(V, ...)
  Z2 <- as.im(mV, na.replace=0)
  imcov(Z, Z2)
}

imcov <- function(X, Y=X) {
  if(missing(Y)) Y <- NULL
  convolve.im(X, Y, reflectX = FALSE, reflectY=TRUE)
}

convolve.im <- function(X, Y=X, ..., reflectX=FALSE, reflectY=FALSE) {
  stopifnot(is.im(X))
  have.Y <- !missing(Y) && !is.null(Y)
  crosscov <- have.Y || reflectX || reflectY
  trap.extra.arguments(..., .Context="In convolve.im")
  #' Check whether Fastest Fourier Transform in the West is available
  west <- fftwAvailable()
  #'
  if(have.Y) {
    # cross-covariance 
    stopifnot(is.im(Y))
    Xbox <- as.rectangle(X)
    Ybox <- as.rectangle(Y)
    # first shift images to common midpoint, to reduce storage
    Xmid <- centroid.owin(Xbox)
    Ymid <- centroid.owin(Ybox)
    svec <- as.numeric(Xmid) - as.numeric(Ymid)
    Y <- shift(Y, svec)
    # ensure images are compatible
    XY <- harmonise.im(X=X, Y=Y)
    X <- XY$X
    Y <- XY$Y
  } else {
    # Y is missing or NULL
    Y <- X
    Xbox <- Ybox <- as.rectangle(X)
  }
  M <- X$v
  M[is.na(M)] <- 0
  xstep <- X$xstep
  ystep <- X$ystep
  # pad with zeroes
  nr <- nrow(M)
  nc <- ncol(M)
  Mpad <- matrix(0, ncol=2*nc, nrow=2*nr)
  Mpad[1:nr, 1:nc] <- M
  lengthMpad <- 4 * nc * nr
  fM <- fft2D(Mpad, west=west)
  if(!crosscov) {
    # compute convolution square
    G <- fft2D(fM^2, inverse=TRUE, west=west)/lengthMpad
  } else {
    # compute set cross-covariance or convolution by FFT
    N <- Y$v
    N[is.na(N)] <- 0
    Npad <- matrix(0, ncol=2*nc, nrow=2*nr)
    Npad[1:nr, 1:nc] <- N
    fN <- fft2D(Npad, west=west)
    if(reflectY) fN <- Conj(fN)
    if(reflectX) fM <- Conj(fM)
    G <- fft2D(fM * fN, inverse=TRUE, west=west)/lengthMpad
  }
#  cat(paste("maximum imaginary part=", max(Im(G)), "\n"))
  G <- Mod(G) * xstep * ystep
  if(reflectX != reflectY) {
    # Currently G[i,j] corresponds to a vector shift of
    #     dy = (i-1) mod nr, dx = (j-1) mod nc.
    # Rearrange this periodic function so that 
    # the origin of translations (0,0) is at matrix position (nr,nc)
    # NB this introduces an extra row and column
    G <- G[ ((-nr):nr) %% (2 * nr) + 1, (-nc):nc %% (2*nc) + 1]
  }
  # Determine spatial domain of full raster image
  XB <- as.rectangle(X)
  YB <- as.rectangle(Y)
  # undo shift
  if(have.Y) YB <- shift(YB, -svec)
  # reflect
  if(reflectX) XB <- reflect(XB)
  if(reflectY) YB <- reflect(YB)
  # Minkowski sum of covering boxes
  xran <- XB$xrange + YB$xrange
  yran <- XB$yrange + YB$yrange
  # Declare spatial domain
  out <- im(G, xrange = xran, yrange=yran)
  if(crosscov) {
    # restrict to actual spatial domain of function
    if(reflectX) Xbox <- reflect(Xbox)
    if(reflectY) Ybox <- reflect(Ybox)
   # Minkowski sum 
    xran <- Xbox$xrange + Ybox$xrange
    yran <- Xbox$yrange + Ybox$yrange   
    XYbox <- owin(xran, yran)
    out <- out[XYbox, rescue=TRUE]
  }
  return(out)
}

