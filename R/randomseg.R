#
# randomseg.R
#
# $Revision: 1.12 $ $Date: 2016/12/01 09:32:41 $
#

rpoisline <- function(lambda, win=owin()) {
  win <- as.owin(win)
  # determine circumcircle
  xr <- win$xrange
  yr <- win$yrange
  xmid <- mean(xr)
  ymid <- mean(yr)
  width <- diff(xr)
  height <- diff(yr)
  rmax <- sqrt(width^2 + height^2)/2
  boundbox <- owin(xmid + c(-1,1) * rmax, ymid + c(-1,1) * rmax)
  # generate poisson lines through circumcircle
  n <- rpois(1, lambda * 2 * pi * rmax)
  if(n == 0) {
    X <- psp(numeric(0), numeric(0), numeric(0), numeric(0),
             marks=integer(0), 
             window=win)
    attr(X, "lines") <- infline(p=numeric(0), theta=numeric(0))
    attr(X, "linemap") <- integer(0)
    return(X)
  }
  theta <- runif(n, max= 2 * pi)
  p <- runif(n, max=rmax)
  # compute intersection points with circle
  q <- sqrt(rmax^2 - p^2)
  co <- cos(theta)
  si <- sin(theta)
  X <- psp(x0= xmid + p * co + q * si,
           y0= ymid + p * si - q * co,
           x1= xmid + p * co - q * si,
           y1= ymid + p * si + q * co,
           marks = seq_len(n),
           window=boundbox, check=FALSE)
  # infinite lines
  L <- infline(p = p + xmid * co + ymid * si,
               theta = theta)
  # clip to window
  X <- X[win]
  # append info
  linemap <- as.integer(marks(X))
  X <- unmark(X)
  attr(X, "lines") <- L
  attr(X, "linemap") <- linemap
  return(X)
}

rlinegrid <- function(angle=45, spacing=0.1, win=owin()) {
  win <- as.owin(win)
  # determine circumcircle
  width <- diff(win$xrange)
  height <- diff(win$yrange)
  rmax <- sqrt(width^2 + height^2)/2
  xmid <- mean(win$xrange)
  ymid <- mean(win$yrange)
  # generate randomly-displaced grid of lines through circumcircle
  u <- runif(1, min=0, max=spacing) - rmax
  if(u >= rmax)   
    return(psp(numeric(0), numeric(0), numeric(0), numeric(0),
               window=win, check=FALSE))
  p <- seq(from=u, to=rmax, by=spacing)
  # compute intersection points with circle
  q <- sqrt(rmax^2 - p^2)
  theta <- pi * ((angle - 90)/180)
  co <- cos(theta)
  si <- sin(theta)
  X <- psp(x0= xmid + p * co + q * si,
           y0= ymid + p * si - q * co,
           x1= xmid + p * co - q * si,
           y1= ymid + p * si + q * co,
           window=owin(xmid+c(-1,1)*rmax, ymid+c(-1,1)*rmax), check=FALSE)
  # clip to window
  X <- X[win]
  return(X)
}
