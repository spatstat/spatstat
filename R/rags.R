#'
#'      rags.R
#'
#'   Alternating Gibbs Sampler
#'
#'      $Revision: 1.4 $  $Date: 2016/10/17 07:18:58 $
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
                          types=NULL, bmax=NULL, ncycles=100) {
  # validate beta by generating first proposal points
  X <- rmpoispp(lambda=beta, lmax=bmax, ..., types=types)
  ntypes <- length(levels(marks(X)))
  check.nmatrix(hradii, ntypes, things="types of points")
  if(any(is.finite(dh <- diag(hradii)) & dh > 0))
    stop("Interaction between points of the same type is not permitted")
  X <- X[order(coords(X)$x)]
  Y <- split(X)
  for(cycle in 1:ncycles) {
    for(i in 1:ntypes) {
      Xi <- Y[[i]]
      ok <- TRUE
      for(j in (1:ntypes)[-i]) {
        if(!any(ok)) break;
        ok <- ok & !is.close(Xi, hradii[i,j], Y[[j]], sorted=TRUE)
      }
      Y[[i]] <- Xi[ok]
    }
    X <- rmpoispp(lambda=beta, lmax=bmax, ..., types=types)
  }
  Z <- do.call(superimpose, Y)
  return(Z)
}
