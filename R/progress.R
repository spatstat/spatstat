#
#   progress.R
#
#   $Revision: 1.11 $  $Date: 2015/10/04 05:13:00 $
#
#   progress plots (envelope representations)
#

dclf.progress <- function(X, ...)
  mctest.progress(X, ..., exponent=2)

mad.progress <- function(X, ...)
  mctest.progress(X, ..., exponent=Inf)

mctest.progress <- local({

  ordstat <- function(z, k) { sort(z, decreasing=TRUE, na.last=TRUE)[k] }

  smoothquantile <- function(z, alpha) {
    min(quantile(density(z), 1-alpha), max(z))
  }
  
  silentmax <- function(z) {
    if(all(is.nan(z))) NaN else max(z[is.finite(z)])
  }

  mctest.progress <- function(X, fun=Lest, ...,
                              exponent=1, nrank=1, interpolate=FALSE, alpha) {
    check.1.real(exponent)
    explain.ifnot(exponent >= 0)
    if(missing(fun) && inherits(X, "envelope"))
      fun <- NULL
    Z <- envelopeProgressData(X, fun=fun, ..., exponent=exponent)
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
      apply(devsim, 1, ordstat, k=nrank)
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
    scr <- 1
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
    scr <- sc
  } else stop("Argument scale should be a function")

  if(is.infinite(exponent)) {
    # MAD
    devdata <- cummax(Deviant(obs-reference, alt, sc))
    devsim <- apply(Deviant(sim-reference, alt, sc), 2, cummax)
    if(deflate) {
      devdata <- scr * devdata
      devsim <-  scr * devsim
    }
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
      devdata <- scr * devdata^(1/exponent)
      devsim <-  scr * devsim^(1/exponent)
    }
    testname <- if(exponent == 2) "Diggle-Cressie-Loosmore-Ford test" else
                if(exponent == 1) "Integral absolute deviation test" else
                paste("Integrated", ordinal(exponent), "Power Deviation test")
  }
  result <- list(R=R, devdata=devdata, devsim=devsim, testname=testname,
                 scaleR=scr)
  return(result)
}

dg.progress <- function(X, fun=Lest, ...,   
                        exponent=2, nsim=19, nsimsub=nsim-1, nrank=1, alpha, 
                        interpolate=FALSE) {
  Xname <- short.deparse(substitute(X))
  env.here <- sys.frame(sys.nframe())
  Xismodel <- is.ppm(X) || is.kppm(X) || is.lppm(X) || is.slrm(X)
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
  alphastring <- paste(100 * alpha, "%%", sep="")
  ## generate or extract simulated patterns and functions
  E <- envelope(X, fun=fun, ..., nsim=nsim,
                savepatterns=TRUE, savefuns=TRUE,
                verbose=FALSE,
                envir.simul=env.here)
  ## get progress data
  PD <- envelopeProgressData(E, fun=fun, ..., nsim=nsim,
                             exponent=exponent, 
                             verbose=FALSE)
  ## get first level MC test significance trace
  T1 <- mctest.sigtrace(E, fun=fun, nsim=nsim, 
                        exponent=exponent,
                        interpolate=interpolate,
                        confint=FALSE, verbose=FALSE, ...)
  R    <- T1$R
  phat <- T1$pest
  ## second level traces
  T2list <- lapply(attr(E, "simpatterns"),
                   mctest.sigtrace,
                   fun=fun, nsim=nsimsub, 
                   exponent=exponent,
                   interpolate=interpolate,
                   confint=FALSE, verbose=FALSE, ...)
  phati <- sapply(T2list, getElement, name="pest")
  ## Dao-Genton procedure
  dgcritrank <- 1 + rowSums(phat > phati)
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
            paste(mcname, alphastring, "critical value for %s"),
            "zero"),
          labl=c("R", "%s(R)", "%s[crit](R)", "0"),
          unitname = unitname(X), fname = fname)
  fvnames(p, ".") <- c("obs", "crit")
  fvnames(p, ".s") <- c("zero", "crit")
  return(p)
}

