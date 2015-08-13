#
#  sigtrace.R
#
#  $Revision: 1.4 $  $Date: 2015/05/05 04:20:36 $
#
#  Significance traces 
#

dclf.sigtrace <- function(X, ...) mctest.sigtrace(X, ..., exponent=2)

mad.sigtrace <- function(X, ...) mctest.sigtrace(X, ..., exponent=Inf)

mctest.sigtrace <- function(X, fun=Lest, ..., exponent=1, alpha=0.05) {
  check.1.real(exponent)
  explain.ifnot(exponent >= 0)
  if(missing(fun) && inherits(X, "envelope"))
    fun <- NULL
  Z <- envelopeProgressData(X, fun=fun, ..., exponent=exponent)
  R       <- Z$R
  devdata <- Z$devdata
  devsim  <- Z$devsim
  nsim     <- ncol(devsim)
  datarank <- apply(devdata < devsim, 1, sum) +
              apply(devdata == devsim, 1, sum)/2 + 1
  pvalue <- datarank/(nsim+1)
  # Agresti-Coull confidence interval
  successes <- datarank - 1
  trials    <- nsim
  z <- qnorm(1 - (1-0.95)/2)
  nplus <- trials + z^2
  pplus <- (successes + z^2/2)/nplus
  sigmaplus <- sqrt(pplus * (1-pplus)/nplus)
  lo <- pplus - z * sigmaplus
  hi <- pplus + z * sigmaplus
  # create fv object
  p <- fv(data.frame(R=R, pest=pvalue, alpha=alpha, lo=lo, hi=hi),
          argu="R", ylab = quote(p(R)), valu="pest", fmla = . ~ R, 
          desc = c("Interval endpoint R",
            "calculated p-value %s",
            "threshold for significance",
            "lower 95%% limit for p-value",
            "upper 95%% limit for p-value"),
          labl=c("R", "%s(R)", paste(alpha), "lo(R)", "hi(R)"),
          unitname = unitname(X), fname = "p")
  fvnames(p, ".") <- c("pest", "alpha")
  try(fvnames(p, ".s") <- c("lo", "hi"), silent=TRUE)
  return(p)
}

