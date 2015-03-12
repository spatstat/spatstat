#
# quadratresample.R
#
# resample a point pattern by resampling quadrats
#
# $Revision: 1.5 $  $Date: 2010/11/25 02:58:49 $
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
  V <- lapply(B, function(z) { w <- z$window;
                               c(w$xrange[1], w$yrange[1]) })
  out <- list()
  if(verbose)
    cat("Generating resampled patterns...")
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
      progressreport(i, nsamples)
  }
  if(nsamples == 1)
    return(out[[1]])
  return(as.solist(out))
}

