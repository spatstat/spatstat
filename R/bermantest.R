#
# bermantest.R
#
# Test statistics from Berman (1986)
#
#  $Revision: 1.16 $  $Date: 2014/06/22 03:07:20 $
#
#

# --------- outdated --------

bermantest <- function(...) {
  message("bermantest is out of date; use berman.test")
#  .Deprecated("berman.test", package="spatstat")
  berman.test(...)
}

bermantest.ppp <- function(...) {
    message("bermantest.ppp is out of date; use berman.test.ppp")
#  .Deprecated("berman.test.ppp", package="spatstat")
  berman.test.ppp(...)
}

bermantest.ppm <- function(...) {
    message("bermantest.ppm is out of date; use berman.test.ppm")
#  .Deprecated("berman.test.ppm", package="spatstat")
  berman.test.ppm(...)
}

bermantest.lpp <- function(...) {
    message("bermantest.lpp is out of date; use berman.test.lpp")
#  .Deprecated("berman.test.lpp", package="spatstat")
  berman.test.lpp(...)
}

bermantest.lppm <- function(...) {
    message("bermantest.lppm is out of date; use berman.test.lppm")
#  .Deprecated("berman.test.lppm", package="spatstat")
  berman.test.lppm(...)
}

# ---------------------------

berman.test <- function(...) {
  UseMethod("berman.test")
}

berman.test.ppp <-
  function(X, covariate,
           which=c("Z1", "Z2"),
           alternative=c("two.sided", "less", "greater"),
           ...) {
    Xname <- short.deparse(substitute(X))
    covname <- short.deparse(substitute(covariate))
    if(is.character(covariate)) covname <- covariate
    which <- match.arg(which)
    alternative <- match.arg(alternative)

    do.call("bermantestEngine",
            resolve.defaults(list(ppm(X), covariate, which, alternative),
                             list(...),
                             list(modelname="CSR",
                                  covname=covname, dataname=Xname)))
}

berman.test.ppm <- function(model, covariate,
                           which=c("Z1", "Z2"),
                           alternative=c("two.sided", "less", "greater"),
                           ...) {
  modelname <- short.deparse(substitute(model))
  covname <- short.deparse(substitute(covariate))
  if(is.character(covariate)) covname <- covariate
  verifyclass(model, "ppm")
  which <- match.arg(which)
  alternative <- match.arg(alternative)
  if(is.poisson(model) && is.stationary(model))
    modelname <- "CSR"
  do.call("bermantestEngine",
          resolve.defaults(list(model, covariate, which, alternative),
                           list(...),
                           list(modelname=modelname,
                                covname=covname,
                                dataname=model$Qname)))
}

berman.test.lpp <-
  function(X, covariate,
           which=c("Z1", "Z2"),
           alternative=c("two.sided", "less", "greater"),
           ...) {
    Xname <- short.deparse(substitute(X))
    covname <- short.deparse(substitute(covariate))
    if(is.character(covariate)) covname <- covariate
    which <- match.arg(which)
    alternative <- match.arg(alternative)

    do.call("bermantestEngine",
            resolve.defaults(list(lppm(X), covariate, which, alternative),
                             list(...),
                             list(modelname="CSR",
                                  covname=covname, dataname=Xname)))
}

berman.test.lppm <- function(model, covariate,
                           which=c("Z1", "Z2"),
                           alternative=c("two.sided", "less", "greater"),
                           ...) {
  modelname <- short.deparse(substitute(model))
  covname <- short.deparse(substitute(covariate))
  if(is.character(covariate)) covname <- covariate
  verifyclass(model, "lppm")
  which <- match.arg(which)
  alternative <- match.arg(alternative)
  if(is.poisson(model) && is.stationary(model))
    modelname <- "CSR"
  do.call("bermantestEngine",
          resolve.defaults(list(model, covariate, which, alternative),
                           list(...),
                           list(modelname=modelname,
                                covname=covname,
                                dataname=model$Xname)))
}

bermantestEngine <- function(model, covariate,
                             which=c("Z1", "Z2"),
                             alternative=c("two.sided", "less", "greater"),
                             ...,
                             modelname, covname, dataname="") {

  csr <- is.poisson(model) && is.stationary(model)
  if(missing(modelname))
    modelname <- if(csr) "CSR" else short.deparse(substitute(model))
  if(missing(covname)) {
    covname <- short.deparse(substitute(covariate))
    if(is.character(covariate)) covname <- covariate
  }

  which <- match.arg(which)
  alternative <- match.arg(alternative)

  if(!is.poisson(model))
    stop("Only implemented for Poisson point process models")

  # ........... first assemble data ...............
  fram <- spatialCDFframe(model, covariate, ...,
                        modelname=modelname,
                        covname=covname,
                        dataname=dataname)
  fvalues <- fram$values
  info    <- fram$info
  # values of covariate at data points
  ZX <- fvalues$ZX
  # transformed to Unif[0,1] under H0
  U  <- fvalues$U
  # values of covariate at pixels
  Zvalues <- fvalues$Zvalues
  # corresponding pixel areas/weights
  weights <- fvalues$weights
  # intensity of model
  lambda  <- fvalues$lambda

  switch(which,
         Z1={
           #......... Berman Z1 statistic .....................
           method <-
             paste("Berman Z1 test of",
                   if(info$csr) "CSR" else "inhomogeneous Poisson process",
                   "in", info$spacename)
           # sum of covariate values at data points
           Sn <- sum(ZX)
           # predicted mean and variance
           lamwt <- lambda * weights
           En    <- sum(lamwt)
           ESn   <- sum(lamwt * Zvalues)
           varSn <- sum(lamwt * Zvalues^2)
           # working, for plot method
           working <- list(meanZX=mean(ZX),
                           meanZ=ESn/En)
           # standardise
           statistic <- (Sn - ESn)/sqrt(varSn)
           names(statistic) <- "Z1"
           p.value <- switch(alternative,
                            two.sided=2 * pnorm(-abs(statistic)),
                            less=pnorm(statistic),
                            greater=pnorm(statistic, lower.tail=FALSE))
           altblurb <- switch(alternative,
                              two.sided="two-sided",
                              less="mean value of covariate at random points is less than predicted under model",
                              greater="mean value of covariate at random points is greater than predicted under model")
           valuename <- paste("covariate",
                              sQuote(paste(covname, collapse="")),
                              "evaluated at points of",
                              sQuote(dataname))
         },
         Z2={
           #......... Berman Z2 statistic .....................
           method <-
             paste("Berman Z2 test of",
                   if(info$csr) "CSR" else "inhomogeneous Poisson process",
                   "in", info$spacename)
           npts <- length(ZX)
           statistic <- sqrt(12/npts) * (sum(U) - npts/2)
           working <- list(meanU=mean(U))
           names(statistic) <- "Z2"
           p.value <- switch(alternative,
                            two.sided=2 * pnorm(-abs(statistic)),
                            less=pnorm(statistic),
                            greater=pnorm(statistic, lower.tail=FALSE))
           altblurb <- switch(alternative,
                              two.sided="two-sided",
                              less="covariate values at random points have lower quantiles than predicted under model",
                              greater="covariate values at random points have higher quantiles than predicted under model")
           valuename <- paste("covariate",
                              sQuote(paste(covname, collapse="")),
                              "evaluated at points of",
                              sQuote(dataname), "\n\t",
                              "and transformed to uniform distribution under",
                              if(info$csr) modelname else sQuote(modelname))
         })
           
  out <- list(statistic=statistic,
              p.value=p.value,
              alternative=altblurb,
              method=method,
              which=which,
              working=working,
              data.name=valuename,
              fram=fram)
  class(out) <- c("htest", "bermantest")
  return(out)
}

plot.bermantest <-
  function(x, ...,
           lwd=par("lwd"), col=par("col"), lty=par("lty"),
           lwd0=lwd, col0=2, lty0=2)
{
  fram <- x$fram
  if(!is.null(fram)) {
    values <- fram$values
    info <- fram$info
  } else {
    # old style
    ks <- x$ks
    values <- attr(ks, "prep")
    info <- attr(ks, "info")
  }
  work <- x$working
  op <- options(useFancyQuotes=FALSE)
  switch(x$which,
         Z1={
           # plot cdf's of Z
           FZ <- values$FZ
           xxx <- get("x", environment(FZ))
           yyy <- get("y", environment(FZ))
           main <- c(x$method,
                     paste("based on distribution of covariate",
                           sQuote(info$covname)),
                     paste("Z1 statistic =", signif(x$statistic, 4)),
                     paste("p-value=", signif(x$p.value, 4)))
           do.call("plot.default",
                   resolve.defaults(
                                    list(x=xxx, y=yyy, type="l"),
                                    list(...),
                                    list(lwd=lwd0, col=col0, lty=lty0),
                                    list(xlab=info$covname,
                                         ylab="probability",
                                         main=main)))
           FZX <- values$FZX
           if(is.null(FZX))
             FZX <- ecdf(values$ZX)
           plot(FZX, add=TRUE, do.points=FALSE, lwd=lwd, col=col, lty=lty)
           abline(v=work$meanZ, lwd=lwd0,col=col0, lty=lty0, xpd=FALSE)
           abline(v=work$meanZX, lwd=lwd,col=col, lty=lty, xpd=FALSE)
         },
         Z2={
           # plot cdf of U
           U <- values$U
           cdfU <- ecdf(U)
           main <- c(x$method,
                     paste("based on distribution of covariate",
                           sQuote(info$covname)),
                     paste("Z2 statistic =", signif(x$statistic, 4)),
                     paste("p-value=", signif(x$p.value, 4)))
           do.call("plot.ecdf",
                   resolve.defaults(
                                    list(cdfU),
                                    list(...),
                                    list(do.points=FALSE, asp=1),
                                    list(lwd=lwd, col=col, lty=lty),
                                    list(xlab="U", ylab="relative frequency"),
                                    list(main=main)))
           abline(0,1,lwd=lwd0,col=col0,lty=lty0, xpd=FALSE)
           abline(v=0.5, lwd=lwd0,col=col0,lty=lty0, xpd=FALSE)
           abline(v=work$meanU, lwd=lwd,col=col,lty=lty, xpd=FALSE)
         })
  options(op)
  return(invisible(NULL))
}



