#'
#'  adaptive.density.R
#'
#'  $Revision: 1.13 $   $Date: 2019/02/04 06:55:53 $
#'

adaptive.density <- function(X, ...) {
  UseMethod("adaptive.density")
}

adaptive.density.ppp <- function(X, f=0.1, ...,
                                 method=c("tile","count"),
                                 fixed=FALSE,
                                 nrep=1, verbose=TRUE) {
  stopifnot(is.ppp(X))
  nX <- npoints(X)
  method <- match.arg(method)
  check.1.real(f)
  if(badprobability(f))
    stop("f should be a probability between 0 and 1")
  check.1.integer(nrep)
  stopifnot(nrep >= 1)
  ##
  ntess <- floor(f * nX)
  if(ntess == 0) {
    ## naive estimate of intensity
    if(f > 0 && verbose)
      splat("Tiny threshold: returning uniform intensity estimate")
    W <- X$window
    lam <- nX/area(W)
    return(as.im(lam, W, ...))
  }
  if(ntess == nX) {
    ## Voronoi/Dirichlet estimate
    tes <- dirichlet(X)
    tesim <- nnmap(X, what="which", ...)
    lam <- 1/tile.areas(tes)
    out <- eval.im(lam[tesim])
    return(out)
  }
  if(nrep > 1) {
    ## estimate is the average of nrep randomised estimates
    total <- 0
    if(verbose)
      cat(paste("Computing", nrep, "intensity estimates..."))
    state <- list()
    for(i in seq_len(nrep)) {
      estimate <- adaptive.density(X, f, ..., method=method, nrep=1)
      total <- eval.im(total + estimate)
      if(verbose) state <- progressreport(i, nrep, state=state)
    }
    if(verbose) cat("Done.\n")
    average <- eval.im(total/nrep)
    return(average)
  }
  ## perform thinning
  if(!fixed) {
    itess <- thinjump(nX, f)
    sampfrac <- f
  } else {
    itess <- sample(seq_len(nX), ntess, replace=FALSE)
    sampfrac <- as.numeric(nX - ntess)/nX
  }
  Xtess <- X[itess]
  ## make tessellation
  tes <- dirichlet(Xtess)
  ## estimate intensity in each tile
  tilemass <- switch(method,
                     tile = 1,
                     count = table(marks(cut(X[-itess], tes))))
  lam <- as.numeric(tilemass)/(sampfrac * tile.areas(tes))
  ## estimate of intensity at each location
  tesim <- nnmap(Xtess, what="which", ...)
  out <- eval.im(lam[tesim])
  return(out)
}
