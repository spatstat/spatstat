#
#  idw.R
#
#  Inverse-distance weighted smoothing
#
#  $Revision: 1.6 $ $Date: 2015/10/21 09:06:57 $

idw <- function(X, power=2, at="pixels", ...) {
  stopifnot(is.ppp(X) && is.marked(X))
  marx <- marks(X)
  if(is.data.frame(marx)) {
    if(ncol(marx) > 1) {
      # multiple columns of marks - process one-by-one
      out <- list()
      for(j in 1:ncol(marx)) 
        out[[j]] <- idw(X %mark% marx[,j], power=power, at=at, ...)
      names(out) <- names(marx)
      switch(at,
             pixels = { out <- as.solist(out) },
             points = { out <- as.data.frame(out) })
      return(out)
    } else 
      marx <- marx[,1]
  }
  if(!is.numeric(marx))
    stop("Marks must be numeric")
  check.1.real(power)
  switch(at,
         pixels = {
           # create grid
           W <- as.mask(as.owin(X), ...)
           dim <- W$dim
           npixels <- prod(dim)
           # call C
           z <- .C("Cidw",
                   x = as.double(X$x),
                   y = as.double(X$y),
                   v = as.double(marx),
                   n = as.integer(npoints(X)),
                   xstart = as.double(W$xcol[1]),
                   xstep  = as.double(W$xstep),
                   nx     = as.integer(dim[2]),
                   ystart = as.double(W$yrow[1]),
                   ystep  = as.double(W$ystep),
                   ny     = as.integer(dim[1]),
                   power  = as.double(power),
                   num    = as.double(numeric(npixels)),
                   den    = as.double(numeric(npixels)),
                   rat    = as.double(numeric(npixels)))
           out <- as.im(matrix(z$rat, dim[1], dim[2]), W=W)
         },
         points={
           npts <- npoints(X)
           z <- .C("idwloo",
                   x = as.double(X$x),
                   y = as.double(X$y),
                   v = as.double(marx),
                   n = as.integer(npts),
                   power  = as.double(power),
                   num    = as.double(numeric(npts)),
                   den    = as.double(numeric(npts)),
                   rat    = as.double(numeric(npts)))
           out <- z$rat
         })
  return(out)
}
