#
#  adaptive.density.R
#
#  $Revision: 1.4 $   $Date: 2011/05/18 01:24:50 $
#
#

adaptive.density <- function(X, f=0.1, ..., nrep=1, verbose=TRUE) {
  stopifnot(is.ppp(X))
  npts <- npoints(X)
  stopifnot(is.numeric(f) && length(f) == 1 && f > 0 & f < 1)
  ntess <- floor(f * npts)
  if(ntess == 0) {
    # naive estimate of intensity
    if(verbose) cat("Tiny threshold: returning uniform intensity estimate")
    W <- X$window
    lam <- npts/area.owin(W)
    return(as.im(lam, W, ...))
  }
  if(nrep > 1) {
    # estimate is the average of nrep randomised estimates
    total <- 0
    if(verbose) cat(paste("Computing", nrep, "intensity estimates..."))
    for(i in seq_len(nrep)) {
      estimate <- adaptive.density(X, f, ..., nrep=1)
      total <- eval.im(total + estimate)
      if(verbose) progressreport(i, nrep)
    }
    if(verbose) cat("Done.\n")
    average <- eval.im(total/nrep)
    return(average)
  }
  ncount <- npts - ntess
  fcount <- ncount/npts
  itess <- sample(seq_len(npts), ntess, replace=FALSE)
  Xtess <- X[itess]
  Xcount <- X[-itess]
  tes <- dirichlet(Xtess)
  meanintensity <- function(x) { x$n/area.owin(x$window) }
  lam <- unlist(lapply(split(Xcount, tes), meanintensity))
  tesim <- as.im(tes, ...)
  out <- eval.im(lam[as.integer(tesim)]/fcount)
  return(out)
}
