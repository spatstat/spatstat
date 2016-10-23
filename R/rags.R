#'
#'      rags.R
#'
#'   Alternating Gibbs Sampler
#'
#'      $Revision: 1.5 $  $Date: 2016/10/23 10:06:38 $
#'
#' Initial implementation for multitype hard core process
#' without interaction within types

rags <- function(model, ..., ncycles=100) {
  if(!is.list(model)) stop("Argument 'model' should be a list")
  if(!all(c("beta", "hradii") %in% names(model)))
    stop("Argument 'model' should have entries 'beta' and 'hradii'")
  do.call(ragsMultiHard, append(model, list(..., ncycles=ncycles)))
}

ragsMultiHard <- function(beta, hradii, ...,
                          types=NULL, bmax=NULL,
                          periodic=FALSE, ncycles=100) {
  ## validate beta by generating first proposal points
  Xprop <- rmpoispp(lambda=beta, lmax=bmax, ..., types=types)
  ntypes <- length(levels(marks(Xprop)))
  check.nmatrix(hradii, ntypes, things="types of points")
  if(any(is.finite(dh <- diag(hradii)) & dh > 0))
    stop("Interaction between points of the same type is not permitted")
  ## initial state empty
  X <- Xprop[integer(0)]
  Y <- split(X)
  ##
  for(cycle in 1:ncycles) {
    if(cycle > 1)
      Xprop <- rmpoispp(lambda=beta, lmax=bmax, ..., types=types)
    Xprop <- Xprop[order(coords(Xprop)$x)]
    Yprop <- split(Xprop)
    for(i in 1:ntypes) {
      Xi <- Yprop[[i]]
      ok <- TRUE
      for(j in (1:ntypes)[-i]) {
        if(!any(ok)) break;
        ok <- ok & !is.close(Xi, hradii[i,j], Y[[j]], sorted=TRUE,
                             periodic=periodic)
      }
      Y[[i]] <- Xi[ok]
    }
  }
  Z <- do.call(superimpose, Y)
  return(Z)
}
