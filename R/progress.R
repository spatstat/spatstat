#
#   progress.R
#
#   $Revision: 1.7 $  $Date: 2014/10/12 00:16:58 $
#
#   progress plots (envelope representations)
#

dclf.progress <- function(X, ..., nrank=1)
  mctest.progress(X, ..., exponent=2, nrank=nrank)

mad.progress <- function(X, ..., nrank=1)
  mctest.progress(X, ..., exponent=Inf, nrank=nrank)

mctest.progress <- local({

  ordstat <- function(z, k) { sort(z, decreasing=TRUE, na.last=TRUE)[k] }
  
  silentmax <- function(z) {
    if(all(is.nan(z))) NaN else max(z[is.finite(z)])
  }

  mctest.progress <- function(X, fun=Lest, ..., exponent=1, nrank=1) {
    check.1.real(exponent)
    explain.ifnot(exponent >= 0)
    if((nrank %% 1) != 0)
      stop("nrank must be an integer")
    if(missing(fun) && inherits(X, "envelope"))
      fun <- NULL
    Z <- envelopeProgressData(X, fun=fun, ..., exponent=exponent)
    R       <- Z$R
    devdata <- Z$devdata
    devsim  <- Z$devsim
    nsim    <- ncol(devsim)
    critval <- if(nrank == 1) apply(devsim, 1, silentmax) else
    apply(devsim, 1, ordstat, k=nrank)
    alpha   <- nrank/(nsim + 1)
    alphastring <- paste(100 * alpha, "%%", sep="")
    # create fv object
    fname  <- if(is.infinite(exponent)) "mad" else
              if(exponent == 2) "T" else paste("D[",exponent,"]", sep="")
    ylab <- if(is.infinite(exponent)) quote(mad(R)) else
            if(exponent == 2) quote(T(R)) else
            eval(substitute(quote(D[p](R)), list(p=exponent)))
    df <- data.frame(R=R, obs=devdata, crit=critval, zero=0)
    p <- fv(df,
            argu="R", ylab=ylab, valu="obs", fmla = . ~ R, 
            desc = c("Interval endpoint R",
              "observed value of test statistic %s",
              paste("Monte Carlo", alphastring, "critical value for %s"),
              "zero"),
            labl=c("R", "%s(R)", "%s[crit](R)", "0"),
            unitname = unitname(X), fname = fname)
    fvnames(p, ".") <- c("obs", "crit")
    fvnames(p, ".s") <- c("zero", "crit")
    return(p)
  }

  mctest.progress
})


# Do not call this function.
# Performs underlying computations

envelopeProgressData <- function(X, fun=Lest, ..., exponent=1,
                                normalize=FALSE, deflate=FALSE) {
  # compute or extract simulated functions
  X <- envelope(X, fun=fun, ..., savefuns=TRUE)
  Y <- attr(X, "simfuns")
  # extract values
  R   <- with(X, .x)
  obs <- with(X, .y)
  reference <- if("theo" %in% names(X)) with(X, theo) else with(X, mmean)
  sim <- as.matrix(as.data.frame(Y))[, -1]
  nsim <- ncol(sim)

  if(is.infinite(exponent)) {
    # MAD
    devdata <- cummax(abs(obs-reference))
    devsim <- apply(abs(sim-reference), 2, cummax)
    testname <- "Maximum absolute deviation test"
  } else {
    dR <- c(0, diff(R))
    a <- (nsim/(nsim - 1))^exponent
    devdata <- a * cumsum(dR * abs(obs - reference)^exponent)
    devsim <- a * apply(dR * abs(sim - reference)^exponent, 2, cumsum)
    if(normalize) {
      devdata <- devdata/R
      devsim <- sweep(devsim, 1, R, "/")
    }
    if(deflate) {
      devdata <- devdata^(1/exponent)
      devsim <- devsim^(1/exponent)
    }
    testname <- if(exponent == 2) "Diggle-Cressie-Loosmore-Ford test" else
                if(exponent == 1) "Integral absolute deviation test" else
                paste("Integrated", ordinal(exponent), "Power Deviation test")
  }
  result <- list(R=R, devdata=devdata, devsim=devsim, testname=testname)
  return(result)
}
