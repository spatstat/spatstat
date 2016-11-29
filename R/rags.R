#'
#'      rags.R
#'
#'   Alternating Gibbs Sampler
#'
#'      $Revision: 1.6 $  $Date: 2016/11/29 05:01:51 $
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
        ok <- ok & !has.close(Xi, hradii[i,j], Y[[j]], sorted=TRUE,
                             periodic=periodic)
      }
      Y[[i]] <- Xi[ok]
    }
  }
  Z <- do.call(superimpose, Y)
  return(Z)
}

ragsAreaInter <- function(beta, eta, r, ...,
                          win=NULL, bmax=NULL,
                          periodic=FALSE, ncycles=100) {
  check.1.real(eta)
  check.1.real(r)
  if(r == 0 || eta == 1) return(rpoispp(beta, win=win, lmax=bmax, ...))
  if(eta < 1)
    stop("Alternating Gibbs algorithm requires eta >= 1", call.=FALSE)
  if(is.function(beta)) {
    beta <- as.im(beta, W=win, ...)
  } else if(is.numeric(beta)) {
    check.1.real(beta)
    stopifnot(beta >= 0)
  } else if(!is.im(beta)) {
    stop("beta should be a number, a pixel image, or a function(x,y)",
         call.=FALSE)
  }
  if(is.im(beta) && is.null(win))
    win <- as.owin(beta)
  kappa <- beta * eta
  loggamma <- log(eta)/(pi * r^2)
  bmax <- if(is.null(bmax)) NULL else c(max(kappa), loggamma)
  B <- if(is.numeric(beta)) c(kappa, loggamma) else
       solist(kappa, as.im(loggamma, W=win))
  H <- matrix(c(0,r,r,0), 2, 2)
  Y <- ragsMultiHard(B, H, types=1:2, bmax=bmax, periodic=periodic,
                     ncycles=ncycles)
  X <- split(Y)[[1]]
  return(X)
}
