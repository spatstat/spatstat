#
# disc.R
#
# $Revision: 1.5 $ $Date: 2013/05/01 05:46:37 $
#
#

disc <- function(radius=1, centre=c(0,0), ..., mask=FALSE, npoly=128) {
  stopifnot(length(centre) == 2)
  stopifnot(length(radius) == 1)
  stopifnot(radius > 0)
  stopifnot(length(npoly) == 1)
  stopifnot(npoly > 2)
  if(!mask) {
    theta <- seq(from=0, to=2*pi, length.out=npoly+1)[-(npoly+1)]
    x <- centre[1] + radius * cos(theta)
    y <- centre[2] + radius * sin(theta)
    W <- owin(poly=list(x=x, y=y))
  } else {
    B <- owin(c(-1,1),c(-1,1))
    B <- as.mask(B, ...)
    indic <- function(x,y,x0,y0,r) as.integer((x-x0)^2 + (y-y0)^2 < r^2)
    IW <- as.im(indic, B, x0=centre[1], y0=centre[2], r=radius)
    W <- levelset(IW, 1, "==")
  }
  return(W)
}

  
