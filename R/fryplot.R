#
#  fryplot.R
#
#  $Revision: 1.6 $ $Date: 2013/04/25 06:37:43 $
#

fryplot <- function(X, ..., width=NULL, from=NULL, to=NULL) {
  Xname <- short.deparse(substitute(X))
  X <- as.ppp(X)
  n <- npoints(X)
  ismarked <- is.marked(X)
  seqn <- seq_len(n)
  from <- if(is.null(from)) seqn else seqn[from]
  to   <- if(is.null(to))   seqn else seqn[to]
  b <- as.rectangle(X)
  halfspan <- with(b, c(diff(xrange), diff(yrange)))/2
  if(!is.null(width)) {
    halfwidth <- ensure2vector(width)/2
    halfspan <- pmin(halfspan, halfwidth)
  }
  bb <- owin(c(-1,1) * halfspan[1], c(-1,1) * halfspan[2])
  do.call("plot.owin",
          resolve.defaults(list(bb),
                           list(...),
                           list(invert=TRUE),
                           list(main=paste("Fry plot of", Xname))))
  xx <- X$x[to]
  yy <- X$y[to]
  if(ismarked) {
    marx <- as.data.frame(marks(X))
    marx <- marx[to, ,drop=FALSE]
  }
  for(i in from) {
    noti <- (to != i)
    dxi <- xx[noti] - xx[i]
    dyi <- yy[noti] - yy[i]
    oki <- (abs(dxi) < halfspan[1]) & (abs(dyi) < halfspan[2])
    if(any(oki)) {
      mki <- if(ismarked) marx[noti, , drop=FALSE] else NULL
      dXi <- ppp(x=dxi[oki], y=dyi[oki], window=bb,
                 marks=mki[oki,],
                 check=FALSE)
      plot(dXi, add=TRUE, ...)
    }
  }
  return(invisible(NULL))
}

frypoints <- function(X) {
  X <- as.ppp(X)
  b <- as.rectangle(X)
  bb <- owin(c(-1,1) * diff(b$xrange), c(-1,1) * diff(b$yrange))
  n <- X$n
  xx <- X$x
  yy <- X$y
  dx <- outer(xx, xx, "-")
  dy <- outer(yy, yy, "-")
  nondiag <- matrix(TRUE, n, n)
  diag(nondiag) <- FALSE
  DX <- as.vector(dx[nondiag])
  DY <- as.vector(dy[nondiag])
  Fry <- ppp(DX, DY, window=bb, check=FALSE)
  if(is.marked(X)) {
    marx <- as.data.frame(marks(X))
    rowind <- row(nondiag)[nondiag]
    marks(Fry) <- marx[rowind, ]
  }
  return(Fry)
}
