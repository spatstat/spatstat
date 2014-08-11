#
#   deltametric.R
#
#   Delta metric
#
#   $Revision: 1.3 $  $Date: 2010/07/10 11:42:08 $
#

deltametric <- function(A, B, p=2, c=Inf, ...) {
  stopifnot(is.numeric(p) && length(p) == 1 && p > 0)
  # ensure frames are identical
  bb <- bounding.box(as.rectangle(A), as.rectangle(B))
  # enforce identical frames
  A <- rebound(A, bb)
  B <- rebound(B, bb)
  # compute distance functions
  dA <- distmap(A, ...)
  dB <- distmap(B, ...)
  if(!is.infinite(c)) {
    dA <- eval.im(pmin(dA, c))
    dB <- eval.im(pmin(dB, c))
  }
  if(is.infinite(p)) {
    # L^infinity
    Z <- eval.im(abs(dA-dB))
    delta <- summary(Z)$max
  } else {
    # L^p
    Z <- eval.im(abs(dA-dB)^p)
    iZ <- summary(Z)$mean
    delta <- iZ^(1/p)
  }
  return(delta)
}





