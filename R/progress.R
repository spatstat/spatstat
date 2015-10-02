#
#   progress.R
#
#   $Revision: 1.8 $  $Date: 2015/10/02 05:09:42 $
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
                                 alternative=c("two.sided", "less", "greater"),
                                 scale=NULL, 
                                 normalize=FALSE, deflate=FALSE) {
  alt <- alternative <- match.arg(alternative)
  # compute or extract simulated functions
  X <- envelope(X, fun=fun, ..., alternative=alternative, savefuns=TRUE)
  Y <- attr(X, "simfuns")
  # extract values
  R   <- with(X, .x)
  obs <- with(X, .y)
  reference <- if("theo" %in% names(X)) with(X, theo) else with(X, mmean)
  sim <- as.matrix(as.data.frame(Y))[, -1]
  nsim <- ncol(sim)

  ## determine rescaling if any
  if(is.null(scale)) {
    sc <- NULL
  } else if(is.function(scale)) {
    sc <- scale(R)
    sname <- "scale(r)"
    ans <- check.nvector(sc, length(R), things="values of r",
                         fatal=FALSE, vname=sname)
    if(!ans)
      stop(attr(ans, "whinge"), call.=FALSE)
    if(any(bad <- (sc <= 0))) {
      ## issue a warning unless this only happens at r=0
      if(any(bad[R > 0]))
        warning(paste("Some values of", sname, "were negative or zero:",
                      "scale was reset to 1 for these values"),
                call.=FALSE)
      sc[bad] <- 1
    }
  } else stop("Argument scale should be a function")

  if(is.infinite(exponent)) {
    # MAD
    devdata <- cummax(Deviant(obs-reference, alt, sc))
    devsim <- apply(Deviant(sim-reference, alt, sc), 2, cummax)
    testname <- "Maximum absolute deviation test"
  } else {
    dR <- c(0, diff(R))
    a <- (nsim/(nsim - 1))^exponent
    devdata <- a * cumsum(dR * Deviant(obs - reference, alt, sc)^exponent)
    devsim <- a * apply(dR * Deviant(sim - reference, alt, sc)^exponent, 2, cumsum)
    if(normalize) {
      devdata <- devdata/R
      devsim <- sweep(devsim, 1, R, "/")
    }
    if(deflate) {
      scr <- sc %orifnull% 1
      devdata <- scr * devdata^(1/exponent)
      devsim <-  scr * devsim^(1/exponent)
    }
    testname <- if(exponent == 2) "Diggle-Cressie-Loosmore-Ford test" else
                if(exponent == 1) "Integral absolute deviation test" else
                paste("Integrated", ordinal(exponent), "Power Deviation test")
  }
  result <- list(R=R, devdata=devdata, devsim=devsim, testname=testname)
  return(result)
}
