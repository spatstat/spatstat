#
#   bw.ppl.R
#
#   Likelihood cross-validation for kernel smoother of point pattern
#
#   $Revision: 1.2 $ $Date: 2013/08/26 02:34:00 $
#

bw.ppl <- function(X, ..., srange=NULL, ns=32) {
  stopifnot(is.ppp(X))
  if(!is.null(srange)) check.range(srange) else
      srange <- c(min(nndist(X)), diameter(as.owin(X))/2)
  sigma <- exp(seq(log(srange[1]), log(srange[2]), length=ns))
  cv <- numeric(ns)
  for(i in 1:ns) {
    si <- sigma[i]
    lamx <- density(X, sigma=si, at="points", leaveoneout=TRUE)
    lam <- density(X, sigma=si)
    cv[i] <- sum(log(lamx)) - integral.im(lam)
  }
  result <- bw.optim(cv, sigma, iopt=which.max(cv), 
                     creator="bw.ppl",
                     criterion="Likelihood Cross-Validation")
  return(result)
}
