#
#  clickjoin.R
#
# interactive addition/deletion of segments between vertices
#

clickjoin <- function(X, ..., add=TRUE, m=NULL, join=TRUE) {
  verifyclass(X, "ppp")
  if(!(is.logical(join) && length(join) == 1))
    stop("join should be a single logical value")
  plot(X, add=add, pch=16)
  if(is.null(m)) {
    m <- matrix(FALSE, npoints(X), npoints(X))
  } else {
    stopifnot(is.matrix(m) && is.logical(m))
    stopifnot(all(dim(m) == npoints(X)))
    from <- as.vector(row(m)[m])
    to   <- as.vector(col(m)[m])
    with(X, segments(x[from], y[from], x[to], y[to]))
  }
  while(TRUE) {
    twoid <- identify(X, plot=FALSE, n=2)
    n <- length(twoid)
    if(n == 0) break
    if(n == 2) {
      m[twoid[1L],twoid[2L]] <- m[twoid[2L],twoid[1L]] <- join
      lines(X$x[twoid], X$y[twoid], ...)
    }
  }
  return(m)
}
