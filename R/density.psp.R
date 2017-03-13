#
#
#  density.psp.R
#
#  $Revision: 1.8 $    $Date: 2016/12/02 10:41:25 $
#
#

density.psp <- function(x, sigma, ..., edge=TRUE,
                        method=c("FFT", "C", "interpreted")) {
  verifyclass(x, "psp")
  method <- match.arg(method)
  w <- x$window
  n <- x$n
  if(missing(sigma))
    sigma <- 0.1 * diameter(w)
  w <- as.mask(w, ...)
  len <- lengths.psp(x)
  if(n == 0 || all(len == 0))
    return(as.im(0, w))
  #
  ang <- angles.psp(x, directed=TRUE)
  xy <- rasterxy.mask(w)
  xx <- xy$x
  yy <- xy$y
  switch(method,
         interpreted = {
           #' compute matrix contribution from each segment 
           coz <- cos(ang)
           zin <- sin(ang)
           for(i in seq_len(n)) {
             en <- x$ends[i,]
             dx <- xx - en$x0
             dy <- yy - en$y0
             u1 <- dx * coz[i] + dy * zin[i]
             u2 <- - dx * zin[i] + dy * coz[i]
             value <- dnorm(u2, sd=sigma) *
               (pnorm(u1, sd=sigma) - pnorm(u1-len[i], sd=sigma))
             totvalue <- if(i == 1L) value else (value + totvalue)
           }
           dens <- im(totvalue, w$xcol, w$yrow)
         },
         C = {
           #' C implementation of the above
           xs <- x$ends$x0
           ys <- x$ends$y0
           xp <- as.numeric(as.vector(xx))
           yp <- as.numeric(as.vector(yy))
           np <- length(xp)
           z <- .C("segdens",
                   sigma = as.double(sigma),
                   ns = as.integer(n),
                   xs = as.double(xs),
                   ys = as.double(ys),
                   alps = as.double(ang),
                   lens = as.double(len),
                   np = as.integer(np),
                   xp = as.double(xp), 
                   yp = as.double(yp),
                   z = as.double(numeric(np)),
                   PACKAGE = "spatstat")
           dens <- im(z$z, w$xcol, w$yrow)
         },
         FFT = {
           L <- pixellate(x, ...)
           L <- L/with(L, xstep * ystep)
           dens <- blur(L, sigma, normalise=edge, bleed=FALSE)
         })
  unitname(dens) <- unitname(x)
  if(edge && method != "FFT") {
    edg <- second.moment.calc(midpoints.psp(x), sigma, what="edge", ...)
    dens <- eval.im(dens/edg)
  }
  dens <- dens[x$window, drop=FALSE]
  return(dens)
}
