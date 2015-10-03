#
#  sigtrace.R
#
#  $Revision: 1.5 $  $Date: 2015/10/03 10:41:01 $
#
#  Significance traces 
#

dclf.sigtrace <- function(X, ...) mctest.sigtrace(X, ..., exponent=2)

mad.sigtrace <- function(X, ...) mctest.sigtrace(X, ..., exponent=Inf)

mctest.sigtrace <- local({
  
  mctest.sigtrace <- function(X, fun=Lest, ..., exponent=1,
                              interpolate=FALSE, alpha=0.05,
                              confint=TRUE) {
    check.1.real(exponent)
    explain.ifnot(exponent >= 0)
    if(missing(fun) && inherits(X, "envelope"))
      fun <- NULL
    Z <- envelopeProgressData(X, fun=fun, ..., exponent=exponent)
    R       <- Z$R
    devdata <- Z$devdata
    devsim  <- Z$devsim
    nsim     <- ncol(devsim)
    if(!interpolate) {
      #' Monte Carlo p-value
      datarank <- apply(devdata < devsim, 1, sum) +
        apply(devdata == devsim, 1, sum)/2 + 1
      pvalue <- datarank/(nsim+1)
    } else {
      #' interpolated p-value
      devs <- cbind(devdata, devsim)
      pvalue <- apply(devs, 1, rowwise.interp.tailprob)
    }
    if(!confint) {
      #' create fv object without confidence interval
      p <- fv(data.frame(R=R, pest=pvalue, alpha=alpha), 
              argu="R", ylab = quote(p(R)), valu="pest", fmla = . ~ R, 
              desc = c("Interval endpoint R",
                "calculated p-value %s",
                "threshold for significance"), 
              labl=c("R", "%s(R)", paste(alpha)), 
              unitname = unitname(X), fname = "p")
      fvnames(p, ".") <- c("pest", "alpha")
    } else {
      # confidence interval
      if(!interpolate) {
        #' Agresti-Coull confidence interval
        successes <- datarank - 1
        trials    <- nsim
        z <- qnorm(1 - (1-0.95)/2)
        nplus <- trials + z^2
        pplus <- (successes + z^2/2)/nplus
        sigmaplus <- sqrt(pplus * (1-pplus)/nplus)
        lo <- pplus - z * sigmaplus
        hi <- pplus + z * sigmaplus
      } else {
        #' confidence interval by delta method
        pSE    <- apply(devs, 1, rowwise.se)
        z <- qnorm(1 - (1-0.95)/2)
        lo <- pmax(0, pvalue - z * pSE)
        hi <- pmin(1, pvalue + z * pSE)
      }
      #' create fv object with confidence interval
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
      fvnames(p, ".s") <- c("lo", "hi")
    }
    return(p)
  }

  ## interpolated p-value
  interpol.tailprob <- function(x, q) {
    sigma <- bw.nrd0(x)
    mean(pnorm(q, mean=x, sd=sigma, lower.tail=FALSE))
  }
  rowwise.interp.tailprob <- function(x) {
    interpol.tailprob(x[-1], x[1])
  }
  ## estimated SE of p-value
  interpol.se <- function(x, q) {
    sigma <- bw.nrd0(x)
    z <- density(x, sigma)
    v <- mean(z$y * pnorm(q, mean=z$x, sd=sigma, lower.tail=FALSE)^2) * diff(range(z$x))
    sqrt(v)/length(x)
  }
  rowwise.se <- function(x) {
    interpol.se(x[-1], x[1])
  }
    
  mctest.sigtrace
})

  
  
