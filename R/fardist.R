##
##  fardist.R
##
## Farthest distance to boundary
##
##  $Revision: 1.11 $ $Date: 2017/06/05 10:31:58 $

fardist <- function(X, ...) {
  UseMethod("fardist")
}

fardist.owin <- function(X, ..., squared=FALSE) {
  verifyclass(X, "owin")
  M <- as.mask(X, ...)
  V <- if(is.mask(X)) vertices(M) else vertices(X)
  nx <- dim(M)[2L]
  ny <- dim(M)[1L]
  x0 <- M$xcol[1L]
  y0 <- M$yrow[1L]
  xstep <- M$xstep
  ystep <- M$ystep
  if(squared) {
    z <- .C("fardist2grid",
            nx = as.integer(nx),
            x0 = as.double(x0),
            xstep = as.double(xstep),
            ny = as.integer(ny),
            y0 = as.double(y0),
            ystep = as.double(ystep),
            np = as.integer(length(V$x)),
            xp = as.double(V$x),
            yp = as.double(V$y),
            dfar = as.double(numeric(nx * ny)),
            PACKAGE = "spatstat")
  } else {
    z <- .C("fardistgrid",
            nx = as.integer(nx),
            x0 = as.double(x0),
            xstep = as.double(xstep),
            ny = as.integer(ny),
            y0 = as.double(y0),
            ystep = as.double(ystep),
            np = as.integer(length(V$x)),
            xp = as.double(V$x),
            yp = as.double(V$y),
            dfar = as.double(numeric(nx * ny)),
            PACKAGE = "spatstat")
  }
  out <- im(z$dfar, xcol=M$xcol, yrow=M$yrow,
            xrange=M$xrange, yrange=M$yrange, unitname=unitname(M))
  if(!is.rectangle(X))
    out <- out[X, drop=FALSE]
  return(out)
}
  
fardist.ppp <- function(X, ..., squared=FALSE) {
  verifyclass(X, "ppp")
  V <- vertices(Window(X))
  D2 <- crossdist(X$x, X$y, V$x, V$y, squared=TRUE) 
  D2max <- apply(D2, 1L, max)
  if(squared) return(D2max) else return(sqrt(D2max))
}

