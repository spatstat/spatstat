#
#  dclftest.R
#
#  $Revision: 1.30 $  $Date: 2015/10/05 06:29:19 $
#
#  Monte Carlo tests for CSR (etc)
#

clf.test <- function(...) {
 .Deprecated("dclf.test", package="spatstat")
 dclf.test(...)
}

dclf.test <- function(X, ...,
                      alternative=c("two.sided", "less", "greater"),
                      rinterval=NULL, scale=NULL,
                      interpolate=FALSE) {
  Xname <- short.deparse(substitute(X))
  envelopeTest(X, ..., exponent=2, alternative=alternative,
                       rinterval=rinterval,
                       scale=scale, interpolate=interpolate,
                       Xname=Xname)
}

mad.test <- function(X, ...,
                     alternative=c("two.sided", "less", "greater"),
                     rinterval=NULL, scale=NULL,
                     interpolate=FALSE) {
  Xname <- short.deparse(substitute(X))
  envelopeTest(X, ..., exponent=Inf, alternative=alternative,
               rinterval=rinterval,
               scale=scale, interpolate=interpolate,
               Xname=Xname)
}

## measure deviation of summary function
## taking account of alternative hypothesis and possible scaling

Deviant <- local({
  
  plusvalue <- function(x) {
    d <- dim(x)
    y <- pmax(0, x)
    if(!is.null(d)) y <- matrix(y, d[1], d[2])
    return(y)
  }

  minusvalue <- function(x) plusvalue(-x)

  Deviant <- function(x, alternative, scaling=NULL) {
    if(!is.null(scaling)) x <- x/scaling
    switch(alternative,
           two.sided = abs(x),
           less = minusvalue(x),
           greater = plusvalue(x))
  }

  Deviant
})

## workhorse function

envelopeTest <-
  function(X, ...,
           exponent=1,
           alternative=c("two.sided", "less", "greater"),
           rinterval=NULL,
           scale=NULL, 
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
    alt <- alternative <- match.arg(alternative)
    force(save.envelope)
    check.1.real(exponent)
    explain.ifnot(exponent >= 0)
    ## compute or extract simulated functions
    X <- envelope(X, ...,
                  savefuns=TRUE, savepatterns=savepatterns,
                  Yname=Xname, verbose=verbose)
    Y <- attr(X, "simfuns")
    ## extract values
    r   <- with(X, .x)
    obs <- with(X, .y)
    sim <- as.matrix(as.data.frame(Y))[, -1]
    nsim <- ncol(sim)
    nr <- length(r)
    ## choose function as reference
    has.theo <- ("theo" %in% names(X))
    use.theo <- identical(attr(X, "einfo")$use.theory, TRUE)
    if(use.theo && !has.theo)
      warning("No theoretical function available; use.theory ignored")
    if(use.theo && has.theo) {
      reference <- with(X, theo)
      theo.used <- TRUE
    } else {
      ## compute sample mean of simulations *and* observed 
      reference <- apply(cbind(sim, obs), 1, mean, na.rm=TRUE)
      theo.used <- FALSE
    }
    ## determine interval of r values for computation
    rok <- r
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
      rok <- r[ok]
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
      sc <- NULL
    } else if(is.function(scale)) {
      sc <- scale(rok)
      sname <- "scale(r)"
      ans <- check.nvector(sc, nr, things="values of r",
                           fatal=FALSE, vname=sname)
      if(!ans)
        stop(attr(ans, "whinge"), call.=FALSE)
      if(any(bad <- (sc <= 0))) {
        ## issue a warning unless this only happens at r=0
        if(any(bad[rok > 0]))
          warning(paste("Some values of", sname, "were negative or zero:",
                        "scale was reset to 1 for these values"),
                  call.=FALSE)
        sc[bad] <- 1
      }
    } else stop("Argument scale should be a function")

    ## compute test statistic
    if(is.infinite(exponent)) {
      ## MAD
      devdata <- max(Deviant(obs-reference, alt, sc))
      names(devdata) <- "mad"
      devsim <- apply(Deviant(sim-reference, alt, sc), 2, max)
      testname <- "Maximum absolute deviation test"
    } else {
      L <- if(nr > 1) diff(rinterval) else 1
      a <- L * (if(theo.used) 1 else ((nsim+1)/nsim)^exponent)
      if(exponent == 2) {
        ## Cramer-von Mises
        devdata <- a * mean((Deviant(obs - reference, alt, sc))^2)
        names(devdata) <- "u"
        devsim <- a * .colMeans((Deviant(sim - reference, alt, sc))^2, nr, nsim)
        testname <- "Diggle-Cressie-Loosmore-Ford test"
      } else if(exponent == 1) {
        ## integral absolute deviation
        devdata <- a * mean(Deviant(obs - reference, alt, sc))
        names(devdata) <- "L1"
        devsim <- a * .colMeans(Deviant(sim - reference, alt, sc), nr, nsim)
        testname <- "Integral absolute deviation test"
      } else {
        ## general p
        devdata <- a * mean(((Deviant(obs - reference, alt, sc))^exponent))
        names(devdata) <- "Lp"
        devsim <- a * .colMeans(((Deviant(sim - reference, alt, sc))^exponent), nr, nsim)
        testname <- paste("Integrated",
                          ordinal(exponent), "Power Deviation test")
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
                          randomise = sample(1:nties, 1))
        datarank <- datarank + tierank
        if(verbose) message("Ties were encountered")
      }
      pvalue <- datarank/(nsim+1)
      ## bookkeeping
      statistic <- data.frame(devdata, rank=datarank)
      colnames(statistic)[1] <- names(devdata)
    } else {
      ## Dao-Genton style interpolation
      fhat <- density(devsim)
      pvalue <- with(fhat, {
        if(max(x) <= devdata) 0 else
        mean(y[x >= devdata]) * (max(x) - devdata)
      })
      statistic <- data.frame(devdata)
      colnames(statistic)[1] <- names(devdata)
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
    testname <- c(paste(testname, "of", nullmodel),
                  paste(testtype, "based on", nsim,
                        "simulations", e$constraints), 
                  paste("Summary function:", fname),
                  paste("Reference function:",
                        if(theo.used) "theoretical" else "sample mean"),
                  paste("Alternative:", alternative),
                  paste("Interval of distance values:",
                        prange(rinterval), uname),
                  scaleblurb
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
                                   interpolate=interpolate,
                                   scale=scale,
                                   tie.rule=tie.rule,
                                   use.theo=use.theo)
    }
    return(result)
  }



    
   
