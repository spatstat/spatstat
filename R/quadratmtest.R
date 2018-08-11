#
#   method for 'quadrat.test' for class mppm
#
#   $Revision: 1.8 $   $Date: 2015/08/12 07:29:17 $
#
quadrat.test.mppm <- function(X, ...) {
  Xname <- short.deparse(substitute(X))
  if(!is.poisson.mppm(X))
    stop("Model is not a Poisson point process")
  
  subs <- subfits(X)
  tests <- anylapply(subs, quadrat.test.ppm, ..., fitname=Xname)

  df.est <- length(coef(X))
  return(pool.quadrattest(tests, Xname=Xname, df.est=df.est))
}

