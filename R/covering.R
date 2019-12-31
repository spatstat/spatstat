#'
#'   covering.R
#'
#'  $Revision: 1.4 $  $Date: 2019/12/31 03:57:05 $
#'

covering <- function(W, r, ..., giveup=1000) {
  W <- as.owin(W)
  ## compute distance to boundary
  D <- distmap(W, invert=TRUE, ...)
  D <- D[W, drop=FALSE]
  M <- as.owin(D)
  pixstep <- max(M$xstep, M$ystep)
  ## very small distances
  if(r <= pixstep) {
    warning("r is smaller than the pixel resolution: returning pixel centres",
            call.=FALSE)
    xy <- rasterxy.mask(M, drop=TRUE)
    return(ppp(xy$x, xy$y, window=W, check=FALSE))
  }
  ## find the point of W farthest from the boundary
  X <- where.max(D)
  ## build a hexagonal grid through this point
  ruse <- if(is.convex(W)) r else (r * 2/3)
  ruse <- max(pixstep, ruse - pixstep)
  H <- hexgrid(W, ruse, offset=c(X$x, X$y), origin=c(0,0))
  if(npoints(H) == 0) H <- X
  ## this may not suffice if W is irregular
  for(i in 1:giveup) {
    DH <- distmap(H)
    if(max(DH) < ruse && npoints(H) > 0) return(H)
    Hnew <- where.max(DH)
    H <- superimpose(H, Hnew, W=W)
  }
  stop(paste("Failed to converge after adding", giveup, "points"), call.=FALSE)
}

