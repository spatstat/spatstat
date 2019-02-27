#
#  dclftest.R
#
#  $Revision: 1.38 $  $Date: 2019/02/27 08:51:51 $
#
#  Monte Carlo tests for CSR (etc)
#

# clf.test <- function(...) {
#  .Deprecated("dclf.test", package="spatstat")
#  dclf.test(...)
# }

dclf.test <- function(X, ...,
                      alternative=c("two.sided", "less", "greater"),
                      rinterval=NULL, leaveout=1, scale=NULL, clamp=FALSE, 
                      interpolate=FALSE) {
  Xname <- short.deparse(substitute(X))
  envelopeTest(X, ..., exponent=2, alternative=alternative,
                       rinterval=rinterval, leaveout=leaveout,
                       scale=scale, clamp=clamp, interpolate=interpolate,
                       Xname=Xname)
}

mad.test <- function(X, ...,
                     alternative=c("two.sided", "less", "greater"),
                     rinterval=NULL, leaveout=1, scale=NULL, clamp=FALSE,
                     interpolate=FALSE) {
  Xname <- short.deparse(substitute(X))
  envelopeTest(X, ..., exponent=Inf, alternative=alternative,
               rinterval=rinterval, leaveout=leaveout, 
               scale=scale, clamp=clamp, interpolate=interpolate,
               Xname=Xname)
}

## measure deviation of summary function
## leaveout = 0: typically 'ref' is exact theoretical value
##               Compute raw deviation.
## leaveout = 1: 'ref' is mean of simulations *and* observed.
##               Use algebra to compute leave-one-out deviation.
## leaveout = 2: 'ref' is mean of simulations
##               Use algebra to compute leave-two-out deviation.

Deviation <- function(x, ref, leaveout, n, xi=x) {
  if(leaveout == 0) return(x-ref)
  if(leaveout == 1) return((x-ref) * (n+1)/n)
  jackmean <- (n * ref - xi)/(n-1)
  return(x - jackmean)
}

## Evaluate signed or absolute deviation,  
## taking account of alternative hypothesis and possible scaling
## (Large positive values always favorable to alternative)

RelevantDeviation <- local({
  
  positivepart <- function(x) {
    d <- dim(x)
    y <- pmax(0, x)
    if(!is.null(d)) y <- matrix(y, d[1L], d[2L])
    return(y)
  }

  negativepart <- function(x) positivepart(-x)

  RelevantDeviation <- function(x, alternative, clamp=FALSE, scaling=NULL) {
    if(!is.null(scaling)) x <- x/scaling
    switch(alternative,
           two.sided = abs(x),
           less = if(clamp) negativepart(x) else -x,
           greater = if(clamp) positivepart(x) else x)
  }

  RelevantDeviation
})

  
## workhorse function

envelopeTest <-
  function(X, ...,
           exponent=1,
           alternative=c("two.sided", "less", "greater"),
           rinterval=NULL,
           leaveout=1,
           scale=NULL,
           clamp=FALSE,
           tie.rule=c("randomise","mean"),
           interpolate=FALSE,
           save.interpolant = TRUE,
           save.envelope = savefuns || savepatterns,
           savefuns = FALSE, 
           savepatterns = FALSE,
           Xname=NULL,
           verbose=TRUE) {
    if(is.null(Xname)) Xname <- short.deparse(substitute(X))
    tie.rule <- match.arg(tie.rule)
    alternative <- match.arg(alternative)
    if(!(leaveout %in% 0:2))
      stop("Argument leaveout should equal 0, 1 or 2")
    force(save.envelope)
    check.1.real(exponent)
    explain.ifnot(exponent >= 0)
    deviationtype <- switch(alternative,
                            two.sided = "absolute",
                            greater = if(clamp) "positive" else "signed",
                            less = if(clamp) "negative" else "signed")
    deviationblurb <- paste(deviationtype, "deviation")
    ## compute or extract simulated functions
    X <- envelope(X, ...,
                  savefuns=TRUE, savepatterns=savepatterns,
                  Yname=Xname, verbose=verbose)
    Y <- attr(X, "simfuns")
    ## extract values
    r   <- with(X, .x)
    obs <- with(X, .y)
    sim <- as.matrix(as.data.frame(Y))[, -1L]
    nsim <- ncol(sim)
    nr <- length(r)
    ## choose function as reference
    has.theo <- ("theo" %in% names(X))
    use.theo <- identical(attr(X, "einfo")$use.theory, TRUE)
    if(use.theo && !has.theo)
      warning("No theoretical function available; use.theory ignored")
    if(use.theo && has.theo) {
      theo.used <- TRUE
      reference <- with(X, theo)
      leaveout <- 0
    } else {
      theo.used <- FALSE
      if(leaveout == 2) {
        ## use sample mean of simulations only
        reference <- apply(sim, 1L, mean, na.rm=TRUE)
      } else {
        ## use sample mean of simulations *and* observed 
        reference <- apply(cbind(sim, obs), 1L, mean, na.rm=TRUE)
      }
    }
    ## determine interval of r values for computation
    rok <- r
    if(!is.null(rinterval)) {
      check.range(rinterval)
      if(max(r) < rinterval[2L]) {
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
      ok <- (rinterval[1L] <= r & r <= rinterval[2L])
      nr <- sum(ok)
      if(nr == 0) {
        ## rinterval is very short: pick nearest r value
        ok <- which.min(abs(r - mean(rinterval)))
        nr <- 1L
      }
      rok <- r[ok]
      obs <- obs[ok]
      sim <- sim[ok, , drop=FALSE]
      reference <- reference[ok]
    } else {
      rinterval <- range(r)
      bad <- !matrowall(is.finite(as.matrix(X)))
      if(any(bad)) {
        if(bad[1L] && !any(bad[-1L])) {
          ## ditch r = 0
          rinterval <- c(r[2L], max(r))
          if(verbose)
            warning(paste("Some function values were infinite or NaN",
                          "at distance r = 0;",
                          "interval of r values was reset to",
                          prange(rinterval)))
          ok <- (rinterval[1L] <= r & r <= rinterval[2L])
          rok <- r[ok]
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

    ## determine rescaling if any
    if(is.null(scale)) {
      scaling <- NULL
    } else if(is.function(scale)) {
      scaling <- scale(rok)
      sname <- "scale(r)"
      ans <- check.nvector(scaling, nr, things="values of r",
                           fatal=FALSE, vname=sname)
      if(!ans)
        stop(attr(ans, "whinge"), call.=FALSE)
      if(any(bad <- (scaling <= 0))) {
        ## issue a warning unless this only happens at r=0
        if(any(bad[rok > 0]))
          warning(paste("Some values of", sname, "were negative or zero:",
                        "scale was reset to 1 for these values"),
                  call.=FALSE)
        scaling[bad] <- 1
      }
    } else stop("Argument scale should be a function")

    ## compute deviations
    rawdevDat <- Deviation(obs, reference, leaveout, nsim, sim[,1L])
    rawdevSim <- Deviation(sim, reference, leaveout, nsim)
    ## evaluate signed/absolute deviation relevant to alternative
    ddat <- RelevantDeviation(rawdevDat, alternative, clamp, scaling)
    dsim <- RelevantDeviation(rawdevSim, alternative, clamp, scaling)

    ## compute test statistic
    if(is.infinite(exponent)) {
      ## MAD
      devdata <- max(ddat)
      devsim <- apply(dsim, 2, max)
      names(devdata) <- "mad"
      testname <- paste("Maximum", deviationblurb, "test")
      statisticblurb <- paste("Maximum", deviationblurb)
    } else {
      L <- if(nr > 1) diff(rinterval) else 1
      if(exponent == 2) {
        ## Cramer-von Mises
        ddat2 <- if(clamp) ddat^2 else (sign(ddat) * ddat^2)
        dsim2 <- if(clamp) dsim^2 else (sign(dsim) * dsim^2)
        devdata <- L * mean(ddat2)
        devsim  <- L * .colMeans(dsim2, nr, nsim)
        names(devdata) <- "u"
        testname <- "Diggle-Cressie-Loosmore-Ford test"
        statisticblurb <- paste("Integral of squared", deviationblurb)
      } else if(exponent == 1) {
        ## integral absolute deviation
        devdata <- L * mean(ddat)
        devsim  <- L * .colMeans(dsim, nr, nsim)
        names(devdata) <- "L1"
        testname <- paste("Integral", deviationblurb, "test")
        statisticblurb <- paste("Integral of", deviationblurb)
      } else {
        ## general p
        if(clamp) {
          ddatp <- ddat^exponent
          dsimp <- dsim^exponent
        } else {
          ddatp <- sign(ddat) * (abs(ddat)^exponent)
          dsimp <- sign(dsim) * (abs(dsim)^exponent)
        }
        devdata <- L * mean(ddatp)
        devsim  <- L * .colMeans(dsimp, nr, nsim)
        names(devdata) <- "Lp"
        testname <- paste("Integrated",
                          ordinal(exponent), "Power Deviation test")
        statisticblurb <- paste("Integral of",
                                ordinal(exponent), "power of",
                                deviationblurb)
      }
    }
    if(!interpolate) {
      ## standard Monte Carlo test 
      ## compute rank and p-value
      datarank <- sum(devdata < devsim) + 1
      nties <- sum(devdata == devsim)
      if(nties > 0) {
        tierank <- switch(tie.rule,
                          mean = nties/2,
                          randomise = sample(1:nties, 1L))
        datarank <- datarank + tierank
        if(verbose) message("Ties were encountered")
      }
      pvalue <- datarank/(nsim+1)
      ## bookkeeping
      statistic <- data.frame(devdata, rank=datarank)
      colnames(statistic)[1L] <- names(devdata)
    } else {
      ## Dao-Genton style interpolation
      fhat <- density(devsim)
      pvalue <- with(fhat, {
        if(max(x) <= devdata) 0 else
        mean(y[x >= devdata]) * (max(x) - devdata)
      })
      statistic <- data.frame(devdata)
      colnames(statistic)[1L] <- names(devdata)
      nties <- 0
    }
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
             func = "model simulated by evaluating function",
             list = "model simulated by drawing patterns from a list",
             "unrecognised model")
    } else "unrecognised model"
    fname <- deparse(attr(X, "ylab"))
    uname <- with(summary(unitname(X)),
                  if(!vanilla) paste(plural, explain) else NULL)
    testtype <- paste0(if(interpolate) "Interpolated " else NULL,
                       "Monte Carlo")
    scaleblurb <- if(is.null(scale)) NULL else
                  paste("Scale function:", paste(deparse(scale), collapse=" "))
    refblurb <- if(theo.used) "theoretical" else "sample mean"
    leaveblurb <- if(leaveout == 0) paste("observed minus", refblurb) else
                  if(leaveout == 1) "leave-one-out" else "leave-two-out"
    testname <- c(paste(testname, "of", nullmodel),
                  paste(testtype, "test based on", nsim,
                        "simulations", e$constraints), 
                  paste("Summary function:", fname),
                  paste("Reference function:", refblurb),
                  paste("Alternative:", alternative),
                  paste("Interval of distance values:",
                        prange(rinterval), uname),
                  scaleblurb,
                  paste("Test statistic:", statisticblurb),
                  paste("Deviation =", leaveblurb)
                  )
    result <- structure(list(statistic = statistic,
                             p.value = pvalue,
                             method = testname,
                             data.name = e$Yname),
                        class="htest")
    attr(result, "rinterval") <- rinterval
    if(save.interpolant && interpolate)
      attr(result, "density") <- fhat
    if(save.envelope) {
      result <- hasenvelope(result, X)
      attr(result, "statistics") <- list(data=devdata, sim=devsim)
      attr(result, "info") <- list(exponent=exponent,
                                   alternative=alternative,
                                   nties=nties,
                                   leaveout=leaveout,
                                   interpolate=interpolate,
                                   scale=scale, clamp=clamp,
                                   tie.rule=tie.rule,
                                   use.theo=use.theo)
    }
    return(result)
  }



    
   
