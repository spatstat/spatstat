#'
#'  densityVoronoi.R
#'
#'  $Revision: 1.16 $   $Date: 2019/02/06 04:46:29 $
#'

densityVoronoi <- function(X, ...) {
  UseMethod("densityVoronoi")
}

densityVoronoi.ppp <- function(X, f=1, ...,
                               counting=FALSE,
                               fixed=FALSE,
                               nrep=1, verbose=TRUE) {
  stopifnot(is.ppp(X))
  nX <- npoints(X)
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
      estimate <- densityVoronoi(X, f, ...,
                                 counting=counting, fixed=fixed, nrep=1)
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
    tessfrac <- f
  } else {
    itess <- sample(seq_len(nX), ntess, replace=FALSE)
    tessfrac <- as.numeric(ntess)/nX
  }
  Xtess <- X[itess]
  ## make tessellation
  tes <- dirichlet(Xtess)
  ## estimate intensity in each tile
  if(!counting) {
    tilemass <- 1
    expansion <- 1/tessfrac
  } else {
    Xcount <- X[-itess]
    tilemap <- tileindex(Xcount$x, Xcount$y, tes)
    tilemass <- as.numeric(table(tilemap))
    expansion <- 1/(1-tessfrac)
  }
  lam <- expansion * tilemass/tile.areas(tes)
  ## estimate of intensity at each location
  tesim <- nnmap(Xtess, what="which", ...)
  out <- eval.im(lam[tesim])
  return(out)
}
