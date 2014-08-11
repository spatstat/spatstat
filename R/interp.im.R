#
# interp.im.R
#
#  $Revision: 1.2 $  $Date: 2007/05/17 16:41:13 $
#

interp.im <- function(Z, x, y) {
  stopifnot(is.im(Z))
  stopifnot(length(x) == length(y))
  if(!is.null(levels(Z)))
    stop("Interpolation is undefined for factor-valued images")
  ok <- inside.owin(x,y, as.owin(Z))
  # get default lookup values (for boundary cases)
  fallback <- Z[ppp(x[ok], y[ok], window=as.rectangle(Z), check=FALSE)]
  # Transform to grid coordinates
  # so that pixel centres are at integer points,
  # bottom left of image is (0,0)
  xx <- (x[ok] - Z$xcol[1])/Z$xstep
  yy <- (y[ok] - Z$yrow[1])/Z$ystep
  # find grid point to left and below
  # (may transgress boundary)
  xlower <- floor(xx)
  ylower <- floor(yy)
  cc <- as.integer(xlower) + 1
  rr <- as.integer(ylower) + 1
  # determine whether (x,y) is above or below antidiagonal in square
  dx <- xx - xlower
  dy <- yy - ylower
  below <- (dx + dy <= 1)
  # if below, interpolate Z(x,y) = (1-x-y)Z(0,0) + xZ(1,0) + yZ(0,1)
  # if above, interpolate Z(x,y) = (x+y-1)Z(1,1) + (1-x)Z(0,1) + (1-y)Z(1,0)
  V <- Z$v
  lukimyu <- function(ccc, rrr, mat, defaults) {
    dimm <- dim(mat)
    within <- (rrr >= 1 & rrr <= dimm[1] & ccc >= 1 & ccc <= dimm[2])
    result <- defaults
    result[within] <- mat[cbind(rrr[within], ccc[within])]
    result
  }
  values <- ifelse(below,
                   ( (1-dx-dy)*lukimyu(cc,rr,V,fallback)
                   + dx*lukimyu(cc+1,rr,V,fallback)
                   + dy*lukimyu(cc,rr+1,V,fallback)
                    ),
                   ( (dx+dy-1)*lukimyu(cc+1,rr+1,V,fallback)
                   + (1-dx)*lukimyu(cc,rr+1,V,fallback)
                   + (1-dy)*lukimyu(cc+1,rr,V,fallback)
                    ))
  result <- numeric(length(x))
  result[ok] <- values
  result[!ok] <- NA
  return(result)
}
