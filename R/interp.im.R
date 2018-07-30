#
# interp.im.R
#
#  $Revision: 1.6 $  $Date: 2018/07/30 14:29:25 $
#

interp.im <- local({

  lukimyu <- function(ccc, rrr, mat, defaults) {
    dimm <- dim(mat)
    within <- (rrr >= 1 & rrr <= dimm[1L] & ccc >= 1 & ccc <= dimm[2L])
    result <- defaults
    result[within] <- mat[cbind(rrr[within], ccc[within])]
    result
  }

  interp.im <- function(Z, x, y=NULL, bilinear=FALSE) {
    stopifnot(is.im(Z))
    if(!is.null(levels(Z)))
      stop("Interpolation is undefined for factor-valued images")
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    ok <- inside.owin(x,y, as.owin(Z))
    V <- Z$v
    ## get default lookup values (for boundary cases)
    fallback <- Z[ppp(x[ok], y[ok], window=as.rectangle(Z), check=FALSE)]
    ## Transform to grid coordinates
    ## so that pixel centres are at integer points,
    ## bottom left of image is (0,0)
    xx <- (x[ok] - Z$xcol[1L])/Z$xstep
    yy <- (y[ok] - Z$yrow[1L])/Z$ystep
    ## find grid point to left and below
    ## (may transgress boundary)
    xlower <- floor(xx)
    ylower <- floor(yy)
    cc <- as.integer(xlower) + 1L
    rr <- as.integer(ylower) + 1L
    dx <- xx - xlower
    dy <- yy - ylower
    if(bilinear) {
      ## 'orthodox'
      values <- ((1-dx) * (1-dy) * lukimyu(cc,rr,V,fallback)
                 + dx   * (1-dy) * lukimyu(cc+1,rr,V,fallback)
                 + (1-dx) * dy   * lukimyu(cc,rr+1,V,fallback)
                 + dx    * dy    * lukimyu(cc+1,rr+1,V,fallback)
      )
    } else {
      ## original & default
      ## determine whether (x,y) is above or below antidiagonal in square
      below <- (dx + dy <= 1)
      ## if below,interpolate Z(x,y) = (1-x-y)Z(0,0) + xZ(1,0) + yZ(0,1)
      ## if above,interpolate Z(x,y) = (x+y-1)Z(1,1) + (1-x)Z(0,1) + (1-y)Z(1,0)
      values <- ifelse(below,
                      ( (1-dx-dy)*lukimyu(cc,rr,V,fallback)
                       + dx*lukimyu(cc+1,rr,V,fallback)
                       + dy*lukimyu(cc,rr+1,V,fallback)
                       ),
                      ( (dx+dy-1)*lukimyu(cc+1,rr+1,V,fallback)
                       + (1-dx)*lukimyu(cc,rr+1,V,fallback)
                       + (1-dy)*lukimyu(cc+1,rr,V,fallback)
                      ))
    }
    result <- numeric(length(x))
    result[ok] <- values
    result[!ok] <- NA
    return(result)
  }

  interp.im
})
