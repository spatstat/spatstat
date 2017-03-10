#
#
#      distmap.R
#
#      $Revision: 1.22 $     $Date: 2017/02/20 06:32:32 $
#
#
#     Distance transforms
#
#
distmap <- function(X, ...) {
  UseMethod("distmap")
}

distmap.ppp <- function(X, ...) {
  verifyclass(X, "ppp")
  e <- exactdt(X, ...)
  W <- e$w
  uni <- unitname(W)
  dmat <- e$d
  imat <- e$i
  V <- im(dmat, W$xcol, W$yrow, unitname=uni)
  I <- im(imat, W$xcol, W$yrow, unitname=uni)
  if(X$window$type == "rectangle") {
    # distance to frame boundary
    bmat <- e$b
    B <- im(bmat, W$xcol, W$yrow, unitname=uni)
  } else {
    # distance to window boundary, not frame boundary
    bmat <- bdist.pixels(W, style="matrix")
    B <- im(bmat, W$xcol, W$yrow, unitname=uni)
    # clip all to window
    V <- V[W, drop=FALSE]
    I <- I[W, drop=FALSE]
    B <- B[W, drop=FALSE]
  }
  attr(V, "index") <- I
  attr(V, "bdry")  <- B
  return(V)
}

distmap.owin <- function(X, ..., discretise=FALSE, invert=FALSE) {
  verifyclass(X, "owin")
  uni <- unitname(X)
  if(X$type == "rectangle") {
    M <- as.mask(X, ...)
    Bdry <- im(bdist.pixels(M, style="matrix"),
               M$xcol, M$yrow, unitname=uni)
    if(!invert)
      Dist <- as.im(M, value=0)
    else 
      Dist <- Bdry
  } else if(X$type == "polygonal" && !discretise) {
    Edges <- edges(X)
    Dist <- distmap(Edges, ...)
    Bdry <- attr(Dist, "bdry")
    if(!invert) 
      Dist[X] <- 0
    else {
      bb <- as.rectangle(X)
      bigbox <- grow.rectangle(bb, diameter(bb)/4)
      Dist[complement.owin(X, bigbox)] <- 0
    }
  } else {
    X <- as.mask(X, ...)
    if(invert)
      X <- complement.owin(X)
    xc <- X$xcol
    yr <- X$yrow
    nr <- X$dim[1L]
    nc <- X$dim[2L]
# pad out the input image with a margin of width 1 on all sides
    mat <- X$m
    pad <- invert # boundary condition is opposite of value inside W
    mat <- cbind(pad, mat, pad)
    mat <- rbind(pad, mat, pad)
# call C routine
    res <- .C("distmapbin",
              xmin=as.double(X$xrange[1L]),
              ymin=as.double(X$yrange[1L]),
              xmax=as.double(X$xrange[2L]),
              ymax=as.double(X$yrange[2L]),
              nr = as.integer(nr),
              nc = as.integer(nc),
              inp = as.integer(as.logical(t(mat))),
              distances = as.double(matrix(0, ncol = nc + 2, nrow = nr + 2)),
              boundary = as.double(matrix(0, ncol = nc + 2, nrow = nr + 2)),
              PACKAGE = "spatstat")
  # strip off margins again
    dist <- matrix(res$distances,
                   ncol = nc + 2, byrow = TRUE)[2:(nr + 1), 2:(nc +1)]
    bdist <- matrix(res$boundary,
                    ncol = nc + 2, byrow = TRUE)[2:(nr + 1), 2:(nc +1)]
  # cast as image objects
    Dist <- im(dist,  xc, yr, unitname=uni)
    Bdry <- im(bdist, xc, yr, unitname=uni)
  }
  attr(Dist, "bdry")  <- Bdry
  return(Dist)
}

distmap.psp <- function(X, ...) {
  verifyclass(X, "psp")
  W <- as.mask(Window(X), ...)
  uni <- unitname(W)
  rxy <- rasterxy.mask(W)
  xp <- rxy$x
  yp <- rxy$y
  E <- X$ends
  big <- 2 * diameter(Frame(W))^2
  z <- NNdist2segments(xp, yp, E$x0, E$y0, E$x1, E$y1, big)
  xc <- W$xcol
  yr <- W$yrow
  Dist <- im(array(sqrt(z$dist2), dim=W$dim), xc, yr, unitname=uni)
  Indx <- im(array(z$index, dim=W$dim), xc, yr, unitname=uni)
  Bdry <- im(bdist.pixels(W, style="matrix"), xc, yr, unitname=uni)
  attr(Dist, "index") <- Indx
  attr(Dist, "bdry")  <- Bdry
  return(Dist)
}

