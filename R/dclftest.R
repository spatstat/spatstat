#
#  dclftest.R
#
#  $Revision: 1.21 $  $Date: 2014/02/20 11:47:40 $
#
#  Monte Carlo tests for CSR (etc)
#

clf.test <- function(...) {
  .Deprecated("dclf.test", package="spatstat")
  dclf.test(...)
}

dclf.test <- function(X, ..., rinterval=NULL, use.theo=FALSE) {
  Xname <- short.deparse(substitute(X))
  envelopeTest(X, ..., power=2,
               use.theo=use.theo, rinterval=rinterval, Xname=Xname)
}

mad.test <- function(X, ..., rinterval=NULL, use.theo=FALSE) {
  Xname <- short.deparse(substitute(X))
  envelopeTest(X, ..., power=Inf,
               use.theo=use.theo, rinterval=rinterval, Xname=Xname)
}

envelopeTest <- function(X, ...,
                         power=1, rinterval=NULL,
                         use.theo=FALSE,
                         tie.rule=c("randomise","mean"),
                         save.envelope = savefuns || savepatterns,
                         savefuns = FALSE, 
                         savepatterns = FALSE, 
                         Xname=NULL,
                         verbose=TRUE,
                         internal=NULL) {
  if(is.null(Xname)) Xname <- short.deparse(substitute(X))
  tie.rule <- match.arg(tie.rule)
  force(save.envelope)
  check.1.real(power)
  explain.ifnot(power >= 0)
  if(use.theo) {
    # using theoretical function as reference.
    # ensure resulting envelope object includes theoretical function.
    internal <- resolve.defaults(internal, list(csr=TRUE))
  }
  # case where X is a previous result of dclf.test, etc
  if(inherits(X, "htest")) {
    if(is.null(envX <- attr(X, "envelope")))
      stop(paste(Xname, "does not contain simulation data"))
    X <- envX
  }
  # compute or extract simulated functions
  X <- envelope(X, ...,
                savefuns=TRUE, savepatterns=savepatterns,
                Yname=Xname, internal=internal, verbose=verbose)
  Y <- attr(X, "simfuns")
  # extract values
  r   <- with(X, .x)
  obs <- with(X, .y)
  sim <- as.matrix(as.data.frame(Y))[, -1]
  nsim <- ncol(sim)
  nr <- length(r)
  # choose function as reference
  has.theo <- ("theo" %in% names(X))
  if(use.theo && !has.theo)
    warning("No theoretical function available; use.theo ignored")
  if(use.theo && has.theo) {
    reference <- with(X, theo)
    used.theo <- TRUE
  } else {
    # compute sample mean of simulations *and* observed 
    reference <- apply(cbind(sim, obs), 1, mean, na.rm=TRUE)
    used.theo <- FALSE
  }
  # determine interval of r values for computation
  if(!is.null(rinterval)) {
    stopifnot(is.numeric(rinterval))
    stopifnot(length(rinterval) == 2)
    stopifnot(rinterval[1] < rinterval[2])
    if(max(r) < rinterval[2]) {
      oldrinterval <- rinterval
      rinterval <- intersect.ranges(rinterval, range(r))
      if(verbose)
        warning(paste("The interval", prange(oldrinterval),
                      "is too large for the available data;",
                      "it has been trimmed to", prange(rinterval)))
    }
    ok <- (rinterval[1] <= r & r <= rinterval[2])
    obs <- obs[ok]
    sim <- sim[ok, ]
    reference <- reference[ok]
  } else {
    rinterval <- range(r)
    bad <- !apply(is.finite(as.matrix(X)), 1, all)
    if(any(bad)) {
      if(bad[1] && !any(bad[-1])) {
        # ditch r = 0
        rinterval <- c(r[2], max(r))
        if(verbose)
          warning(paste("Some function values were infinite or NaN",
                        "at distance r = 0; interval of r values was reset to",
                        prange(rinterval)))
        ok <- (rinterval[1] <= r & r <= rinterval[2])
        obs <- obs[ok]
        sim <- sim[ok, ]
        reference <- reference[ok]
      } else {
        # problem
        rbadmax <- max(r[bad])
        unitinfo <- summary(unitname(X))
        stop(paste("Some function values were infinite or NaN",
                   "at distances r up to",
                   paste(rbadmax, ".", sep=""),
                   "Please specify a shorter", sQuote("rinterval")))
      }
    } 
  }

  # compute test statistic
  if(is.infinite(power)) {
    # MAD
    devdata <- max(abs(obs-reference))
    names(devdata) <- "mad"
    devsim <- apply(abs(sim-reference), 2, max)
    testname <- "Maximum absolute deviation test"
  } else {
    a <- diff(rinterval) * (if(used.theo) 1 else ((nsim+1)/nsim)^power)
    if(power == 2) {
      # Cramer-von Mises
      devdata <- a * mean((obs - reference)^2)
      names(devdata) <- "u"
      devsim <- a * .colMeans((sim - reference)^2, nr, nsim)
      testname <- "Diggle-Cressie-Loosmore-Ford test"
    } else if(power == 1) {
      # integral absolute deviation
      devdata <- a * mean(abs(obs - reference))
      names(devdata) <- "L1"
      devsim <- a * .colMeans(abs(sim - reference), nr, nsim)
      testname <- "Integral absolute deviation test"
    } else {
      # general p
      devdata <- a * mean((abs(obs - reference)^power))
      names(devdata) <- "Lp"
      devsim <- a * .colMeans((abs(sim - reference)^power), nr, nsim)
      testname <- paste("Integrated", ordinal(power), "Power Deviation test")
    }
  }
  # compute rank and p-value
  datarank <- sum(devdata < devsim) + 1
  nties <- sum(devdata == devsim)
  if(nties > 0) {
    tierank <- switch(tie.rule,
                      mean = nties/2,
                      randomise = sample(1:nties, 1))
    datarank <- datarank + tierank
    if(verbose) message("Ties were encountered")
  }
  pvalue <- datarank/(nsim+1)
  # bookkeeping
  statistic <- data.frame(devdata, rank=datarank)
  colnames(statistic)[1] <- names(devdata)
  e <- attr(X, "einfo")
  nullmodel <-
    if(identical(e$csr, TRUE)) "CSR" else 
    if(!is.null(e$simtype)) {
      switch(e$simtype,
             csr = "CSR",
             rmh = paste("fitted",
               if(identical(e$pois, TRUE)) "Poisson" else "Gibbs",
               "model"),
             kppm = "fitted cluster model",
             expr = "model simulated by evaluating expression",
             list = "model simulated by drawing patterns from a list",
             "unrecognised model")
    } else "unrecognised model"
  fname <- deparse(attr(X, "ylab"))
  uname <- with(summary(unitname(X)),
                if(!vanilla) paste(plural, explain) else NULL)
  testname <- c(paste(testname, "of", nullmodel),
                paste("Monte Carlo test based on", nsim, "simulations"),
                paste("Summary function:", fname),
                paste("Reference function:",
                      if(used.theo) "theoretical" else "sample mean"),
                paste("Interval of distance values:",
                      prange(rinterval), uname)
                )
  result <- structure(list(statistic = statistic,
                           p.value = pvalue,
                           method = testname,
                           data.name = e$Yname),
                      class="htest")
  attr(result, "rinterval") <- rinterval
  if(save.envelope) {
    attr(result, "envelope") <- X
    attr(result, "statistics") <- list(data=devdata, sim=devsim)
    attr(result, "info") <- list(power=power,
                                 nties=nties,
                                 tie.rule=tie.rule,
                                 use.theo=use.theo)
  }
  return(result)
}


    
   
