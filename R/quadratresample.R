#
# quadratresample.R
#
# resample a point pattern by resampling quadrats
#
# $Revision: 1.7 $  $Date: 2015/10/21 09:06:57 $
#

quadratresample <- function(X, nx, ny=nx, ...,
                            replace=FALSE, nsamples=1,
                            verbose=(nsamples > 1)) {
  stopifnot(is.ppp(X))
  if(X$window$type != "rectangle")
    stop("Resampling is only implemented for rectangular windows")
  # create tessellation
  A <- quadrats(X, nx=nx, ny=ny)
  # split data over tessellation
  B <- split(X, A)
  nq <- length(B)
  # determine bottom left corner of each tile
  V <- lapply(B, framebottomleft)
  out <- list()
  if(verbose) {
    cat("Generating resampled patterns...")
    pstate <- list()
  }
  for(i in 1:nsamples) {
    # resample tiles
    ind <- sample(1:nq, nq, replace=replace)
    Xresampled <- X
    Bresampled <- B
    for(j in 1:nq) {
      k <- ind[j]
      Bresampled[[j]] <- shift(B[[k]], unlist(V[[j]]) - unlist(V[[k]]))
    }
    split(Xresampled, A) <- Bresampled
    out[[i]] <- Xresampled
    if(verbose)
      pstate <- progressreport(i, nsamples, state=pstate)
  }
  if(nsamples == 1)
    return(out[[1]])
  return(as.solist(out))
}

