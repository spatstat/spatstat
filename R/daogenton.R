#
#  daogenton.R
#
#  Dao-Genton adjusted p-values
#
#  $Revision: 1.6 $  $Date: 2015/08/22 10:08:30 $
#

dg.test <- function(X, ..., exponent=2, nsim=19, nsimsub=nsim-1,
                    alternative=c("two.sided", "less", "greater"),
                    interpolate=FALSE,
                    verbose=TRUE) {
  Xname <- short.deparse(substitute(X))
  alternative <- match.arg(alternative)
  env.here <- sys.frame(sys.nframe())
  Xismodel <- is.ppm(X) || is.kppm(X) || is.lppm(X) || is.slrm(X)
  if(!missing(nsimsub) && !relatively.prime(nsim, nsimsub))
    stop("nsim and nsimsub must be relatively prime")
  # top-level test
  if(verbose) cat("Applying test to original data... ")
  tX <- envelopeTest(X, ...,
                     nsim=nsim, alternative=alternative,
                     interpolate=interpolate,
                     exponent=exponent, savepatterns=TRUE, verbose=FALSE,
                     envir.simul=env.here)
  if(verbose) cat("Done.\n")
  pX <- tX$p.value
  # extract simulated patterns 
  Ylist <- attr(attr(tX, "envelope"), "simpatterns")
  # apply same test to each simulated pattern
  if(verbose) cat(paste("Running tests on", nsim, "simulated patterns... "))
  pY <- numeric(nsim)
  for(i in 1:nsim) {
    if(verbose) progressreport(i, nsim)
    Yi <- Ylist[[i]]
    # if X is a model, fit it to Yi. Otherwise the implicit model is CSR.
    if(Xismodel) Yi <- update(X, Yi)
    tYi <- envelopeTest(Yi, ...,
                        nsim=nsimsub, alternative=alternative,
                        interpolate=interpolate,
                        exponent=exponent, savepatterns=TRUE, verbose=FALSE,
                        envir.simul=env.here)
    pY[i] <- tYi$p.value
  }
  # compute adjusted p-value
  padj <- (1 + sum(pY <= pX))/(1+nsim)
  # pack up
  method <- tX$method
  method <- c("Dao-Genton adjusted goodness-of-fit test",
              paste("based on", method[1]),
              method[-1])
  names(pX) <- "p0"
  result <- structure(list(statistic = pX,
                           p.value = padj,
                           method = method,
                           data.name = Xname),
                      class="htest") 
  attr(result, "rinterval") <- attr(tX, "rinterval")
  attr(result, "pX") <- pX
  attr(result, "pY") <- sort(pY)
  return(result)
}

dg.envelope <- function(X, ..., nsim=19,
                        nsimsub=nsim-1,
                        nrank=1,
                        alternative=c("two.sided", "less", "greater"),
                        verbose=TRUE) {
  Xname <- short.deparse(substitute(X))
  alternative <- match.arg(alternative)
  env.here <- sys.frame(sys.nframe())
  Xismodel <- is.ppm(X) || is.kppm(X) || is.lppm(X) || is.slrm(X)
  # top-level test
  if(verbose) cat("Applying test to original data... ")
  tX <- envelopeTest(X, ...,
                     alternative=alternative,
                     nsim=nsim, nrank=nrank,
                     exponent=Inf, savepatterns=TRUE, savefuns=TRUE,
                     verbose=FALSE,
                     envir.simul=env.here)
  if(verbose) cat("Done.\n")
  ## extract info
  envX <- attr(tX, "envelope")
  ## extract simulated patterns and functions
  SimFuns <- attr(envX, "simfuns")
  Ylist <- attr(envX, "simpatterns")
  # apply same test to each simulated pattern
  if(verbose) cat(paste("Running tests on", nsim, "simulated patterns... "))
  rankY <- numeric(nsim)
  for(i in 1:nsim) {
    if(verbose) progressreport(i, nsim)
    Yi <- Ylist[[i]]
    # if X is a model, fit it to Yi. Otherwise the implicit model is CSR.
    if(Xismodel) Yi <- update(X, Yi)
    tYi <- envelopeTest(Yi, ...,
                        alternative=alternative,
                        nsim=nsimsub, nrank=nrank,
                        exponent=Inf, savepatterns=TRUE, verbose=FALSE,
                        envir.simul=env.here)
    rankY[i] <- tYi$p.value * (nsimsub + 1)
  }
  ## find critical rank 'l' 
  dg.rank <- sort(rankY, na.last=TRUE)[nrank]
  if(verbose) cat("dg.rank=", dg.rank, fill=TRUE)
  ## extract deviation values from top-level simulation
  simdev <- attr(tX, "statistics")[["sim"]]
  ## find critical deviation
  dg.crit <- sort(simdev, decreasing=TRUE, na.last=TRUE)[dg.rank]
  if(verbose) cat("dg.crit=", dg.crit, fill=TRUE)
  alpha <- nrank/nsim
  alphatext <- paste0(100*alpha, "%%")
  ## make fv object, for now
  refname <- if("theo" %in% names(envX)) "theo" else "mmean"
  fname <- attr(envX, "fname")
  result <- (as.fv(envX))[, c(fvnames(envX, ".x"),
                            fvnames(envX, ".y"),
                            refname)]
  refval <- envX[[refname]]
  ## 
  newdata <- data.frame(hi=refval + dg.crit,
                        lo=refval - dg.crit)
  newlabl <- c(makefvlabel(NULL, NULL, fname, "hi"),
               makefvlabel(NULL, NULL, fname, "lo"))
  newdesc <- c(paste("upper", alphatext, "critical boundary for %s"),
               paste("lower", alphatext, "critical boundary for %s"))
  switch(alternative,
         two.sided = { },
         less = {
           newdata$hi <- Inf
           newlabl[1] <- "infinity"
           newdesc[1] <- "infinite upper limit"
         },
         greater = {
           newdata$lo <- -Inf
           newlabl[2] <- "infinity"
           newdesc[2] <- "infinite lower limit"
         })
  result <- bind.fv(result, newdata, newlabl, newdesc)
  fvnames(result, ".") <- rev(fvnames(result, "."))
  fvnames(result, ".s") <- c("lo", "hi")
  return(result)
}

