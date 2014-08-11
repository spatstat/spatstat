#
#   method for 'quadrat.test' for class mppm
#
#   $Revision: 1.7 $   $Date: 2012/09/06 03:50:17 $
#
quadrat.test.mppm <- function(X, ...) {
  Xname <- short.deparse(substitute(X))
  if(!is.poisson.mppm(X))
    stop("Model is not a Poisson point process")
  
  subs <- subfits(X)
  tests <- lapply(subs, quadrat.test.ppm, ..., fitname=Xname)
  class(tests) <- c("listof", class(tests))

  df.est <- length(coef(X))
  return(pool.quadrattest(tests, Xname=Xname, df.est=df.est))
}

