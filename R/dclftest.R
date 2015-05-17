#
#  dclftest.R
#
#  $Revision: 1.24 $  $Date: 2014/08/20 10:39:43 $
#
#  Monte Carlo tests for CSR (etc)
#

clf.test <- function(...) {
 .Deprecated("dclf.test", package="spatstat")
 dclf.test(...)
}

dclf.test <- function(X, ...,
                      alternative=c("two.sided", "less", "greater"),
                      rinterval=NULL, use.theo=FALSE) {
  Xname <- short.deparse(substitute(X))
  envelopeTest(X, ..., exponent=2, alternative=alternative,
                       use.theo=use.theo, rinterval=rinterval, Xname=Xname)
}

mad.test <- function(X, ...,
                     alternative=c("two.sided", "less", "greater"),
                     rinterval=NULL, use.theo=FALSE) {
  Xname <- short.deparse(substitute(X))
  envelopeTest(X, ..., exponent=Inf, alternative=alternative,
               use.theo=use.theo, rinterval=rinterval, Xname=Xname)
}

envelopeTest <- local({

  plusvalue <- function(x) {
    d <- dim(x)
    y <- pmax(0, x)
    if(!is.null(d)) y <- matrix(y, d[1], d[2])
    return(y)
  }

  envelopeTest <-
    function(X, ...,
             exponent=1,
             alternative=c("two.sided", "less", "greater"),
             rinterval=NULL,
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
      alternative <- match.arg(alternative)
      force(save.envelope)
      check.1.real(exponent)
      explain.ifnot(exponent >= 0)
      if(use.theo) {
        ## using theoretical function as reference.
        ## ensure resulting envelope object includes theoretical function.
        internal <- resolve.defaults(internal, list(csr=TRUE))
      }
      ## case where X is a previous result of dclf.test, etc
      if(inherits(X, "htest")) {
        if(is.null(envX <- attr(X, "envelope")))
          stop(paste(Xname, "does not contain simulation data"))
        X <- envX
      }
      ## compute or extract simulated functions
      X <- envelope(X, ...,
                    savefuns=TRUE, savepatterns=savepatterns,
                    Yname=Xname, internal=internal, verbose=verbose)
      Y <- attr(X, "simfuns")
      ## extract values
      r   <- with(X, .x)
      obs <- with(X, .y)
      sim <- as.matrix(as.data.frame(Y))[, -1]
      nsim <- ncol(sim)
      nr <- length(r)
      ## choose function as reference
      has.theo <- ("theo" %in% names(X))
      if(use.theo && !has.theo)
        warning("No theoretical function available; use.theo ignored")
      if(use.theo && has.theo) {
        reference <- with(X, theo)
        used.theo <- TRUE
      } else {
        ## compute sample mean of simulations *and* observed 
        reference <- apply(cbind(sim, obs), 1, mean, na.rm=TRUE)
        used.theo <- FALSE
      }
      ## determine interval of r values for computation
      if(!is.null(rinterval)) {
        check.range(rinterval)
        if(max(r) < rinterval[2]) {
          oldrinterval <- rinterval
          rinterval <- intersect.ranges(rinterval, range(r), fatal=FALSE)
          if(is.null(rinterval))
            stop(paste("The specified rinterval",
                       prange(oldrinterval),
                       "has empty intersection",
                       "with the range of r values",
                       prange(range(r)), 
                       "computed by the summary function"),
                 call.=FALSE)
          if(verbose)
            warning(paste("The interval", prange(oldrinterval),
                          "is too large for the available data;",
                          "it has been trimmed to", prange(rinterval)))
        }
        ok <- (rinterval[1] <= r & r <= rinterval[2])
        nr <- sum(ok)
        if(nr == 0) {
          ## rinterval is very short: pick nearest r value
          ok <- which.min(abs(r - mean(rinterval)))
          nr <- 1
        }
        obs <- obs[ok]
        sim <- sim[ok, , drop=FALSE]
        reference <- reference[ok]
      } else {
        rinterval <- range(r)
        bad <- !apply(is.finite(as.matrix(X)), 1, all)
        if(any(bad)) {
          if(bad[1] && !any(bad[-1])) {
            ## ditch r = 0
            rinterval <- c(r[2], max(r))
            if(verbose)
              warning(paste("Some function values were infinite or NaN",
                            "at distance r = 0;",
                            "interval of r values was reset to",
                            prange(rinterval)))
            ok <- (rinterval[1] <= r & r <= rinterval[2])
            obs <- obs[ok]
            sim <- sim[ok, ]
            reference <- reference[ok]
            nr <- sum(ok)
          } else {
            ## problem
            rbadmax <- paste(max(r[bad]), summary(unitname(X))$plural)
            stop(paste("Some function values were infinite or NaN",
                       "at distances r up to",
                       paste(rbadmax, ".", sep=""),
                       "Please specify a shorter", sQuote("rinterval")))
          }
        } 
      }

      deviant <- switch(alternative,
                        two.sided = function(x) abs(x),
                        less = function(x) plusvalue(-x),
                        greater = plusvalue)

      ## compute test statistic
      if(is.infinite(exponent)) {
        ## MAD
        devdata <- max(deviant(obs-reference))
        names(devdata) <- "mad"
        devsim <- apply(deviant(sim-reference), 2, max)
        testname <- "Maximum absolute deviation test"
      } else {
        L <- if(nr > 1) diff(rinterval) else 1
        a <- L * (if(used.theo) 1 else ((nsim+1)/nsim)^exponent)
        if(exponent == 2) {
          ## Cramer-von Mises
          devdata <- a * mean((deviant(obs - reference))^2)
          names(devdata) <- "u"
          devsim <- a * .colMeans((deviant(sim - reference))^2, nr, nsim)
          testname <- "Diggle-Cressie-Loosmore-Ford test"
        } else if(exponent == 1) {
          ## integral absolute deviation
          devdata <- a * mean(deviant(obs - reference))
          names(devdata) <- "L1"
          devsim <- a * .colMeans(deviant(sim - reference), nr, nsim)
          testname <- "Integral absolute deviation test"
        } else {
          ## general p
          devdata <- a * mean(((deviant(obs - reference))^exponent))
          names(devdata) <- "Lp"
          devsim <- a * .colMeans(((deviant(sim - reference))^exponent), nr, nsim)
          testname <- paste("Integrated",
                            ordinal(exponent), "Power Deviation test")
        }
      }
      ## compute rank and p-value
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
      ## bookkeeping
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
                    paste("Monte Carlo test based on", nsim,
                          "simulations", e$constraints), 
                    paste("Summary function:", fname),
                    paste("Reference function:",
                          if(used.theo) "theoretical" else "sample mean"),
                    paste("Alternative:", alternative),
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
        attr(result, "info") <- list(exponent=exponent,
                                     alternative=alternative,
                                     nties=nties,
                                     tie.rule=tie.rule,
                                     use.theo=use.theo)
      }
      return(result)
    }

  envelopeTest
})



    
   
