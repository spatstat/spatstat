##
## disc.R
##
##  discs and ellipses
##
## $Revision: 1.11 $ $Date: 2014/03/31 05:57:35 $
##

disc <- local({

  indic <- function(x,y,x0,y0,r) { as.integer((x-x0)^2 + (y-y0)^2 < r^2) }
  
  disc <- function(radius=1, centre=c(0,0), ..., mask=FALSE, npoly=128) {
    stopifnot(length(radius) == 1)
    stopifnot(radius > 0)
    centre <- as2vector(centre)
    stopifnot(length(npoly) == 1)
    stopifnot(npoly > 2)
    if(!mask) {
      theta <- seq(from=0, to=2*pi, length.out=npoly+1)[-(npoly+1)]
      x <- centre[1] + radius * cos(theta)
      y <- centre[2] + radius * sin(theta)
      W <- owin(poly=list(x=x, y=y), check=FALSE)
    } else {
      xr <- centre[1] + radius * c(-1,1)
      yr <- centre[2] + radius * c(-1,1)
      B <- owin(xr,yr)
      IW <- as.im(indic, B, x0=centre[1], y0=centre[2], r=radius, ...)
      W <- levelset(IW, 1, "==")
    }
    return(W)
  }

  disc
})

ellipse <- local({
  
  indic <- function(x,y,x0,y0,a,b,co,si){
    x <- x-x0
    y <- y-y0
    as.integer(((x*co + y*si)/a)^2 + ((-x*si + y*co)/b)^2 < 1)
  }

  ellipse <- function(a, b, centre=c(0,0), phi=0, ...,
                      mask=FALSE, npoly = 128) {
    ## Czechs:
    stopifnot(length(a) == 1)
    stopifnot(a > 0)
    stopifnot(length(b) == 1)
    stopifnot(b > 0)
    centre <- as2vector(centre)
    stopifnot(length(phi) == 1)
    stopifnot(length(npoly) == 1)
    stopifnot(npoly > 2)
    ## Rotator cuff:
    co <- cos(phi)
    si <- sin(phi)
    ## Mask:
    if(mask) {
      ## Thetas maximizing x and y.
      tx <- atan(-b*tan(phi)/a)
      ty <- atan(b/(a*tan(phi)))
      ## Maximal x and y (for centre = c(0,0)).
      xm <- a*co*cos(tx) - b*si*sin(tx)
      ym <- a*si*cos(ty) + b*co*sin(ty)
      ## Range of x and y.
      xr <- xm*c(-1,1)+centre[1]
      yr <- ym*c(-1,1)+centre[2]
      ## Wrecked-angle to contain the mask.
      B  <- as.mask(owin(xr,yr),...)
      ## Build the mask as a level set.
      IW <- as.im(indic, B, x0=centre[1], y0=centre[2], a=a, b=b, co=co, si=si)
      return(levelset(IW, 1, "=="))
    }
    ## Polygonal.
    ## Build "horizontal" ellipse centred at 0:
    theta <- seq(0, 2 * pi, length = npoly+1)[-(npoly+1)]
    xh <-  a * cos(theta)
    yh <-  b * sin(theta)

    ## Rotate through angle phi and shift centre:
    x  <- centre[1] + co*xh - si*yh
    y  <- centre[2] + si*xh + co*yh
    owin(poly=list(x = x, y = y))
  }

  ellipse
})

