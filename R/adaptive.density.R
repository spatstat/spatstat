#
#  adaptive.density.R
#
#  $Revision: 1.7 $   $Date: 2015/05/11 10:40:44 $
#
#

adaptive.density <- function(X, f=0.1, ..., nrep=1, verbose=TRUE) {
  stopifnot(is.ppp(X))
  npts <- npoints(X)
  check.1.real(f)
  if(badprobability(f))
    stop("f should be a probability between 0 and 1")
  ntess <- floor(f * npts)
  if(ntess == 0) {
    # naive estimate of intensity
    if(f > 0 && verbose)
      splat("Tiny threshold: returning uniform intensity estimate")
    W <- X$window
    lam <- npts/area(W)
    return(as.im(lam, W, ...))
  }
  if(ntess == npts) {
    ## Voronoi/Dirichlet estimate
    tes <- dirichlet(X)
#    tesim <- as.im(tes, ...)
    tesim <- nnmap(X, what="which", ...)
    lam <- 1/tile.areas(tes)
    out <- eval.im(lam[tesim])
    return(out)
  }
  if(nrep > 1) {
    # estimate is the average of nrep randomised estimates
    total <- 0
    if(verbose)
      cat(paste("Computing", nrep, "intensity estimates..."))
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
  lam <- unlist(lapply(split(Xcount, tes), intensity))
#  tesim <- as.im(tes, ...)
#  out <- eval.im(lam[as.integer(tesim)]/fcount)
  tesim <- nnmap(Xtess, what="which", ...)
  out <- eval.im(lam[tesim]/fcount)
  return(out)
}
