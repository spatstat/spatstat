#
#   progress.R
#
#   $Revision: 1.20 $  $Date: 2015/12/15 10:40:12 $
#
#   progress plots (envelope representations)
#

dclf.progress <- function(X, ...)
  mctest.progress(X, ..., exponent=2)

mad.progress <- function(X, ...)
  mctest.progress(X, ..., exponent=Inf)

mctest.progress <- local({

  smoothquantile <- function(z, alpha) {
    min(quantile(density(z), 1-alpha), max(z))
  }
  
  silentmax <- function(z) {
    if(all(is.nan(z))) return(NaN)
    z <- z[is.finite(z)]
    if(length(z) == 0) return(NA) else return(max(z))
  }

  mctest.progress <- function(X, fun=Lest, ...,
                              exponent=1, nrank=1, interpolate=FALSE,
                              alpha, rmin=0) {
    check.1.real(exponent)
    explain.ifnot(exponent >= 0)
    if(missing(fun) && inherits(X, "envelope"))
      fun <- NULL
    Z <- envelopeProgressData(X, fun=fun, ..., rmin=rmin, exponent=exponent)
    R       <- Z$R
    devdata <- Z$devdata
    devsim  <- Z$devsim
    nsim    <- ncol(devsim)
    # determine 'alpha' and 'nrank'
    if(missing(alpha)) {
      if((nrank %% 1) != 0)
        stop("nrank must be an integer")
      alpha   <- nrank/(nsim + 1)
    } else {
      check.1.real(alpha)
      stopifnot(alpha > 0 && alpha < 1)
      if(!interpolate) {
        if(!missing(nrank))
          warning("nrank was ignored because alpha was given", call.=FALSE)
        nrank <- alpha * (nsim + 1)
        if(abs(nrank - round(nrank)) > 1e-2)
          stop("alpha should be a multiple of 1/(nsim + 1)", call.=FALSE)
        nrank <- as.integer(round(nrank))
      }
    }
    alphastring <- paste(100 * alpha, "%%", sep="")
    # compute critical values
    critval <-
      if(interpolate) apply(devsim, 1, smoothquantile, alpha=alpha) else
      if(nrank == 1) apply(devsim, 1, silentmax) else
      apply(devsim, 1, orderstats, k=nrank, decreasing=TRUE)
    # create fv object
    fname  <- if(is.infinite(exponent)) "mad" else
              if(exponent == 2) "T" else paste("D[",exponent,"]", sep="")
    ylab <- if(is.infinite(exponent)) quote(mad(R)) else
            if(exponent == 2) quote(T(R)) else
            eval(substitute(quote(D[p](R)), list(p=exponent)))
    df <- data.frame(R=R, obs=devdata, crit=critval, zero=0)
    mcname <- if(interpolate) "interpolated Monte Carlo" else "Monte Carlo"
    p <- fv(df,
            argu="R", ylab=ylab, valu="obs", fmla = . ~ R, 
            desc = c("Interval endpoint R",
              "observed value of test statistic %s",
              paste(mcname, alphastring, "critical value for %s"),
              "zero"),
            labl=c("R", "%s(R)", "%s[crit](R)", "0"),
            unitname = unitname(X), fname = fname)
    fvnames(p, ".") <- c("obs", "crit", "zero")
    fvnames(p, ".s") <- c("zero", "crit")
    p <- hasenvelope(p, Z$envelope)  # envelope may be NULL
    return(p)
  }

  mctest.progress
})


# Do not call this function.
# Performs underlying computations

envelopeProgressData <- local({
  envelopeProgressData <-
    function(X, fun=Lest, ..., exponent=1,
             alternative=c("two.sided", "less", "greater"),
             leaveout=1, scale=NULL, clamp=FALSE, 
             normalize=FALSE, deflate=FALSE,
             rmin=0,
             save.envelope = savefuns || savepatterns,
             savefuns = FALSE, 
             savepatterns = FALSE) {
    alternative <- match.arg(alternative)
    if(!(leaveout %in% 0:2))
      stop("Argument leaveout should equal 0, 1 or 2")
    ## compute or extract simulated functions
    X <- envelope(X, fun=fun, ..., alternative=alternative,
                  savefuns=TRUE, savepatterns=savepatterns)
    Y <- attr(X, "simfuns")
    ## extract values
    R   <- with(X, .x)
    obs <- with(X, .y)
    sim <- as.matrix(as.data.frame(Y))[, -1]
    nsim <- ncol(sim)
    ## choose function as reference
    has.theo <- ("theo" %in% names(X))
    use.theo <- identical(attr(X, "einfo")$use.theory, TRUE)
    if(use.theo && !has.theo)
      warning("No theoretical function available; use.theory ignored")
    if(use.theo && has.theo) {
#      theo.used <- TRUE
      reference <- with(X, theo)
      leaveout <- 0
    } else {
#      theo.used <- FALSE
      if(leaveout == 2) {
        ## use sample mean of simulations only
        reference <- with(X, mmean)
      } else {
        ## use sample mean of simulations *and* observed 
        reference <- (nsim * with(X, mmean) + obs)/(nsim + 1)
      }
    }
    ## restrict range
    if(rmin > 0) {
      if(sum(R >= rmin) < 2)
        stop("rmin is too large for the available range of r values")
      nskip <- sum(R < rmin)
    } else nskip <- 0
  
    ## determine rescaling if any
    if(is.null(scale)) {
      scaling <- NULL
      scr <- 1
    } else if(is.function(scale)) {
      scaling <- scale(R)
      sname <- "scale(r)"
      ans <- check.nvector(scaling, length(R), things="values of r",
                           fatal=FALSE, vname=sname)
      if(!ans)
        stop(attr(ans, "whinge"), call.=FALSE)
      if(any(bad <- (scaling <= 0))) {
        ## issue a warning unless this only happens at r=0
        if(any(bad[R > 0]))
          warning(paste("Some values of", sname, "were negative or zero:",
                        "scale was reset to 1 for these values"),
                  call.=FALSE)
        scaling[bad] <- 1
      }
      scr <- scaling
    } else stop("Argument scale should be a function")

    ## compute deviations
    rawdevDat <- Deviation(obs, reference, leaveout, nsim, sim[,1])
    rawdevSim <- Deviation(sim, reference, leaveout, nsim)
    ## evaluate signed/absolute deviation relevant to alternative
    ddat <- RelevantDeviation(rawdevDat, alternative, clamp, scaling)
    dsim <- RelevantDeviation(rawdevSim, alternative, clamp, scaling)

    ## compute test statistics
    if(is.infinite(exponent)) {
      ## MAD
      devdata <- cummaxskip(ddat, nskip)
      devsim <- apply(dsim, 2, cummaxskip, nskip=nskip)
      if(deflate) {
        devdata <- scr * devdata
        devsim <-  scr * devsim
      }
      testname <- "Maximum absolute deviation test"
    } else {
      dR <- c(0, diff(R))
      if(clamp || (alternative == "two.sided")) {
        ## deviations are nonnegative
        devdata <- cumsumskip(dR * ddat^exponent, nskip)
        devsim  <- apply(dR * dsim^exponent, 2, cumsumskip, nskip=nskip)
      } else {
        ## sign of deviations should be retained
        devdata <- cumsumskip(dR * sign(ddat) * abs(ddat)^exponent,
                                  nskip=nskip)
        devsim  <- apply(dR * sign(dsim) * abs(dsim)^exponent,
                         2, cumsumskip, nskip=nskip)
      }
      if(normalize) {
        devdata <- devdata/R
        devsim <- sweep(devsim, 1, R, "/")
      }
      if(deflate) {
        devdata <- scr * sign(devdata) * abs(devdata)^(1/exponent) 
        devsim <-  scr * sign(devsim) * abs(devsim)^(1/exponent) 
      }
      testname <- if(exponent == 2) "Diggle-Cressie-Loosmore-Ford test" else
                  if(exponent == 1) "Integral absolute deviation test" else
                  paste("Integrated", ordinal(exponent), "Power Deviation test")
    }
    result <- list(R=R, devdata=devdata, devsim=devsim, testname=testname,
                   scaleR=scr, clamp=clamp)
    if(save.envelope) 
      result$envelope <- X
    return(result)
  }

  cumsumskip <- function(x, nskip=0) {
    if(nskip == 0) cumsum(x) else c(rep(NA, nskip), cumsum(x[-seq_len(nskip)]))
  }

  cummaxskip <- function(x, nskip=0) {
    if(nskip == 0) cummax(x) else c(rep(NA, nskip), cummax(x[-seq_len(nskip)]))
  }

  envelopeProgressData
})

dg.progress <- function(X, fun=Lest, ...,   
                        exponent=2, nsim=19, nsimsub=nsim-1, nrank=1, alpha, 
                        leaveout=1, interpolate=FALSE, rmin=0, 
                        savefuns=FALSE, savepatterns=FALSE,
                        verbose=TRUE) {
  env.here <- sys.frame(sys.nframe())
  if(!missing(nsimsub) && !relatively.prime(nsim, nsimsub))
    stop("nsim and nsimsub must be relatively prime")
  ## determine 'alpha' and 'nrank'
  if(missing(alpha)) {
    if((nrank %% 1) != 0)
      stop("nrank must be an integer")
    alpha   <- nrank/(nsim + 1)
  } else {
    check.1.real(alpha)
    stopifnot(alpha > 0 && alpha < 1)
    if(!interpolate) {
      if(!missing(nrank))
        warning("nrank was ignored because alpha was given", call.=FALSE)
      nrank <- alpha * (nsim + 1)
      if(abs(nrank - round(nrank)) > 1e-2)
        stop("alpha should be a multiple of 1/(nsim + 1)", call.=FALSE)
      nrank <- as.integer(round(nrank))
    }
  }
  if(verbose)
    cat("Computing first-level test data...")
  ## generate or extract simulated patterns and functions
  E <- envelope(X, fun=fun, ..., nsim=nsim,
                savepatterns=TRUE, savefuns=TRUE,
                verbose=FALSE,
                envir.simul=env.here)
  ## get progress data
  PD <- envelopeProgressData(E, fun=fun, ..., rmin=rmin, nsim=nsim,
                             exponent=exponent, leaveout=leaveout,
                             verbose=FALSE)
  ## get first level MC test significance trace
  T1 <- mctest.sigtrace(E, fun=fun, nsim=nsim, 
                        exponent=exponent,
                        leaveout=leaveout,
                        interpolate=interpolate, rmin=rmin,
                        confint=FALSE, verbose=FALSE, ...)
  R    <- T1$R
  phat <- T1$pest
  if(verbose) {
    cat("Done.\nComputing second-level data... ")
    state <- list()
  }
  ## second level traces
  simpat <- attr(E, "simpatterns")
  phat2 <- matrix(, length(R), nsim)
  for(j in seq_len(nsim)) {
    simj <- simpat[[j]]
    sigj <- mctest.sigtrace(simj,
                            fun=fun, nsim=nsimsub, 
                            exponent=exponent,
                            interpolate=interpolate,
                            leaveout=leaveout,
                            rmin=rmin,
                            confint=FALSE, verbose=FALSE, ...)
    phat2[,j] <- sigj$pest
    if(verbose) state <- progressreport(j, nsim, state=state)
  }
  if(verbose) cat("Done.\n")
  ## Dao-Genton procedure
  dgcritrank <- 1 + rowSums(phat > phat2)
  dgcritrank <- pmin(dgcritrank, nsim)
  devsim.sort <- t(apply(PD$devsim, 1, sort, decreasing=TRUE, na.last=TRUE))
  ii <- cbind(seq_along(dgcritrank), dgcritrank)
  devcrit <- devsim.sort[ii]
  devdata <- PD$devdata
  ## create fv object
  fname  <- if(is.infinite(exponent)) "mad" else
            if(exponent == 2) "T" else paste("D[",exponent,"]", sep="")
  ylab <- if(is.infinite(exponent)) quote(mad(R)) else
          if(exponent == 2) quote(T(R)) else
          eval(substitute(quote(D[p](R)), list(p=exponent)))
  df <- data.frame(R=R, obs=devdata, crit=devcrit, zero=0)
  mcname <- if(interpolate) "interpolated Monte Carlo" else "Monte Carlo"
  p <- fv(df,
          argu="R", ylab=ylab, valu="obs", fmla = . ~ R, 
          desc = c("Interval endpoint R",
            "observed value of test statistic %s",
            paste(mcname, paste0(100 * alpha, "%%"), "critical value for %s"),
            "zero"),
          labl=c("R", "%s(R)", "%s[crit](R)", "0"),
          unitname = unitname(X), fname = fname)
  fvnames(p, ".") <- c("obs", "crit", "zero")
  fvnames(p, ".s") <- c("zero", "crit")
  if(savefuns || savepatterns)
    p <- hasenvelope(p, E)
  return(p)
}

