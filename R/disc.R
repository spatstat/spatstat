##
## disc.R
##
##  discs and ellipses
##
## $Revision: 1.18 $ $Date: 2017/01/15 05:25:16 $
##

disc <- local({

  indic <- function(x,y,x0,y0,r) { as.integer((x-x0)^2 + (y-y0)^2 < r^2) }
  
  disc <- function(radius=1, centre=c(0,0), ...,
                   mask=FALSE, npoly=128, delta=NULL) {
    check.1.real(radius)
    stopifnot(radius > 0)
    centre <- as2vector(centre)
    if(!missing(npoly) && !is.null(npoly) && !is.null(delta))
      stop("Specify either npoly or delta")
    if(!missing(npoly) && !is.null(npoly)) {
      stopifnot(length(npoly) == 1)
      stopifnot(npoly >= 3)
    } else if(!is.null(delta)) {
      check.1.real(delta)
      stopifnot(delta > 0)
      npoly <- pmax(16, ceiling(2 * pi * radius/delta))
    } else npoly <- 128
    if(!mask) {
      theta <- seq(from=0, to=2*pi, length.out=npoly+1)[-(npoly+1L)]
      x <- centre[1L] + radius * cos(theta)
      y <- centre[2L] + radius * sin(theta)
      W <- owin(poly=list(x=x, y=y), check=FALSE)
    } else {
      xr <- centre[1L] + radius * c(-1,1)
      yr <- centre[2L] + radius * c(-1,1)
      B <- owin(xr,yr)
      IW <- as.im(indic, B, x0=centre[1L], y0=centre[2L], r=radius, ...)
      W <- levelset(IW, 1, "==")
    }
    return(W)
  }

  disc
})

hexagon <- function(edge=1, centre=c(0,0), ...,
                    align=c("bottom", "top", "left", "right", "no")) {
  regularpolygon(6, edge, centre, align=align)
}

regularpolygon <- function(n, edge=1, centre=c(0,0), ...,
                           align=c("bottom", "top", "left", "right", "no")) {
  check.1.integer(n)
  check.1.real(edge)
  stopifnot(n >= 3)
  stopifnot(edge > 0)
  align <- match.arg(align)
  theta <- 2 * pi/n
  radius <- edge/(2 * sin(theta/2))
  result <- disc(radius, centre, npoly=n, mask=FALSE)
  if(align != "no") {
    k <- switch(align,
                bottom = 3/4,
                top = 1/4,
                left = 1/2,
                right = 1)
    alpha <- theta * (1/2 - (k * n) %% 1)
    result <- rotate(result, -alpha)
  }
  Frame(result) <- boundingbox(result)
  return(result)
}


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
      xr <- xm*c(-1,1)+centre[1L]
      yr <- ym*c(-1,1)+centre[2L]
      ## Wrecked-angle to contain the mask.
      B  <- as.mask(owin(xr,yr),...)
      ## Build the mask as a level set.
      IW <- as.im(indic, B, x0=centre[1L], y0=centre[2L], a=a, b=b, co=co, si=si)
      return(levelset(IW, 1, "=="))
    }
    ## Polygonal.
    ## Build "horizontal" ellipse centred at 0:
    theta <- seq(0, 2 * pi, length = npoly+1)[-(npoly+1L)]
    xh <-  a * cos(theta)
    yh <-  b * sin(theta)

    ## Rotate through angle phi and shift centre:
    x  <- centre[1L] + co*xh - si*yh
    y  <- centre[2L] + si*xh + co*yh
    owin(poly=list(x = x, y = y))
  }

  ellipse
})

