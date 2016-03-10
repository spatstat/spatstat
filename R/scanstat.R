##
##  scanstat.R
##
##  Spatial scan statistics
##
##  $Revision: 1.15 $  $Date: 2016/02/11 10:17:12 $
##

scanmeasure <- function(X, ...){
  UseMethod("scanmeasure")
}


scanmeasure.ppp <- function(X, r, ..., method=c("counts", "fft")) {
  method <- match.arg(method)
  check.1.real(r)
  ## enclosing window
  R <- as.rectangle(as.owin(X))
  ## determine pixel resolution  
  M <- as.mask(R, ...)
  ## expand domain to include centres of all circles intersecting R
  W <- grow.mask(M, r)
  ## 
  switch(method,
         counts = {
           ## direct calculation using C code
           ## get new dimensions
           dimyx <- W$dim
           xr <- W$xrange
           yr <- W$yrange
           nr <- dimyx[1]
           nc <- dimyx[2]
           ##
           n <- npoints(X)
           zz <- .C("scantrans",
                    x=as.double(X$x),
                    y=as.double(X$y),
                    n=as.integer(n),
                    xmin=as.double(xr[1]),
                    ymin=as.double(yr[1]),
                    xmax=as.double(xr[2]),
                    ymax=as.double(yr[2]),
                    nr=as.integer(nr),
                    nc=as.integer(nc),
                    R=as.double(r),
                    counts=as.integer(numeric(prod(dimyx))))
           zzz <- matrix(zz$counts, nrow=dimyx[1], ncol=dimyx[2], byrow=TRUE)
           Z <- im(zzz, xrange=xr, yrange=yr, unitname=unitname(X))
         },
         fft = {
           ## Previous version of scanmeasure.ppp had
           ##    Y <- pixellate(X, ..., padzero=TRUE)
           ## but this is liable to Gibbs phenomena.
           ## Instead, convolve with small Gaussian (sd = 1 pixel width)
           sigma <- with(W, unique(c(xstep, ystep)))
           Y <- density(X, ..., sigma=sigma)
           ## invoke scanmeasure.im
           Z <- scanmeasure(Y, r)
           Z <- eval.im(as.integer(round(Z)))
         })
  return(Z)
}

scanmeasure.im <- function(X, r, ...) {
  D <- disc(radius=r)
  eps <- with(X, c(xstep,ystep))
  if(any(eps >= 2 * r)) return(eval.im(X * pi * r^2))
  D <- as.im(as.mask(D, eps=eps))
  Z <- imcov(X, D)
  return(Z)
}

scanPoisLRTS <- function(nZ, nG, muZ, muG, alternative) {
  nZco <- nG - nZ
  muZco <- muG - muZ
  nlogn <- function(n, a) ifelse(n == 0, 0, n * log(n/a))
  ll <- nlogn(nZ, muZ) + nlogn(nZco, muZco) - nlogn(nG, muG)
  criterion <- (nZ * muZco - muZ * nZco)
  switch(alternative,
         less={
           ll[criterion > 0] <- 0
         },
         greater={
           ll[criterion < 0] <- 0
         },
         two.sided={})
  return(2 * ll)
}

scanBinomLRTS <- function(nZ, nG, muZ, muG, alternative) {
  nZco <- nG - nZ
  muZco <- muG - muZ
  nlogn <- function(n, a) ifelse(n == 0, 0, n * log(n/a))
  logbin <- function(k, n) { nlogn(k, n) + nlogn(n-k, n) }
  ll <- logbin(nZ, muZ) + logbin(nZco, muZco) - logbin(nG, muG)
  criterion <- (nZ * muZco - muZ * nZco)
  switch(alternative,
         less={
           ll[criterion > 0] <- 0
         },
         greater={
           ll[criterion < 0] <- 0
         },
         two.sided={})
  return(2 * ll)
}

scanLRTS <- function(X, r, ...,
                       method=c("poisson", "binomial"),
                       baseline=NULL,
                       case=2,
                       alternative=c("greater", "less", "two.sided"),
                       saveopt = FALSE,
                       Xmask=NULL) {
  stopifnot(is.ppp(X))
  stopifnot(check.nvector(r))
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  if(is.null(Xmask)) Xmask <- as.mask(as.owin(X), ...)
  switch(method,
         poisson={
           Y <- X
           if(is.null(baseline)) {
             mu <- as.im(Xmask, value=1)
           } else if(is.ppm(baseline)) {
             if(is.marked(baseline))
               stop("baseline is a marked point process: not supported")
             mu <- predict(baseline, locations=Xmask)
           } else if(is.im(baseline) || is.function(baseline)) {
             mu <- as.im(baseline, W=Xmask)
           } else stop(paste("baseline should be",
                             "a pixel image, a function, or a fitted model"))
           nG <- npoints(Y)
         },
         binomial={
           stopifnot(is.multitype(X))
           lev <- levels(marks(X))
           if(length(lev) != 2)
             warning("X should usually be a bivariate (2-type) point pattern")
           if(!is.null(baseline))
             stop("baseline is not supported in the binomial case")
           if(is.character(case) && !(case %in% lev))
             stop(paste("Unrecognised label for cases:", sQuote(case)))
           if(is.numeric(case) && !(case %in% seq_along(lev)))
             stop(paste("Undefined level:", case))
           Y <- split(X)[[case]]
           nG <- npoints(Y)
           mu <- unmark(X)
         })
  ## The following line ensures that the same pixel resolution information
  ## is passed to the two calls to 'scanmeasure' below
  Y$window <- Xmask
  ## 
  nr <- length(r)
  lrts <- vector(mode="list", length=nr)
  for(i in 1:nr) {
    ri <- r[i]
    nZ <- scanmeasure(Y, ri)
    muZ <- scanmeasure(mu, ri)
    if(!compatible.im(nZ, muZ)) {
      ha <- harmonise.im(nZ, muZ)
      nZ <- ha[[1]]
      muZ <- ha[[2]]
    }
    switch(method,
           poisson = {
             muG <- integral.im(mu)
             lrts[[i]] <- eval.im(scanPoisLRTS(nZ, nG, muZ, muG, alternative))
           },
           binomial = {
             muG <- npoints(mu)
             lrts[[i]] <- eval.im(scanBinomLRTS(nZ, nG, muZ, muG, alternative))
           })
  }
  if(length(lrts) == 1) {
    result <- lrts[[1]]
  } else {
    result <- im.apply(lrts, max)
    if(saveopt)
      attr(result, "iopt") <- im.apply(lrts, which.max)
  }
  return(result)
}

scan.test <- function(X, r, ...,
                      method=c("poisson", "binomial"),
                      nsim = 19,
                      baseline=NULL,
                      case = 2,
                      alternative=c("greater", "less", "two.sided"),
                      verbose=TRUE) {
  dataname <- short.deparse(substitute(X))
  stopifnot(is.ppp(X))
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  stopifnot(is.numeric(r))
  check.1.real(nsim)
  if(!(round(nsim) == nsim && nsim > 1))
    stop("nsim should be an integer > 1")
  regionname <-
    paste("circles of radius",
          if(length(r) == 1) r else paste("between", min(r), "and", max(r)))
  ##
  ## compute observed loglikelihood function
  ## This also validates the arguments.
  obsLRTS <- scanLRTS(X=X, r=r,
                          method=method,
                          alternative=alternative, baseline=baseline,
                          case=case, ..., saveopt=TRUE)
  obs <- max(obsLRTS)
  sim <- numeric(nsim)
  ## determine how to simulate
  switch(method,
         binomial={
           methodname <- c("Spatial scan test",
                           "Null hypothesis: constant relative risk",
                           paste("Candidate cluster regions:", regionname),
                           "Likelihood: binomial",
                           paste("Monte Carlo p-value based on",
                                 nsim, "simulations"))

           lev <- levels(marks(X))
           names(lev) <- lev
           casename <- lev[case]
           counted <- paste("points with mark", sQuote(casename), 
                            "inside cluster region")
           simexpr <- expression(rlabel(X))
         },
         poisson={
           counted <- paste("points inside cluster region")
           X <- unmark(X)
           Xwin <- as.owin(X)
           Xmask <- as.mask(Xwin, ...)
           if(is.null(baseline)) {
             nullname <- "Complete Spatial Randomness (CSR)"
             lambda <- intensity(X)
             simexpr <- expression(runifpoispp(lambda, Xwin))
             dont.complain.about(lambda)
           } else if(is.ppm(baseline)) {
             nullname <- baseline$callstring
             rmhstuff <- rmh(baseline, preponly=TRUE, verbose=FALSE)
             simexpr <- expression(rmhEngine(rmhstuff))
             dont.complain.about(rmhstuff)
           } else if(is.im(baseline) || is.function(baseline)) {
             nullname <- "Poisson process with intensity proportional to baseline"
             base <- as.im(baseline, W=Xmask)
             alpha <- npoints(X)/integral.im(base)
             lambda <- eval.im(alpha * base)
             simexpr <- expression(rpoispp(lambda))
             dont.complain.about(lambda)
           } else stop(paste("baseline should be",
                             "a pixel image, a function, or a fitted model"))
           methodname <- c("Spatial scan test",
                           paste("Null hypothesis:", nullname),
                           paste("Candidate cluster regions:", regionname),
                           "Likelihood: Poisson",
                           paste("Monte Carlo p-value based on",
                                 nsim, "simulations"))
         })
  if(verbose) {
    cat("Simulating...")
    pstate <- list()
  }
  for(i in 1:nsim) {
    if(verbose) pstate <- progressreport(i, nsim, state=pstate)
    Xsim <- eval(simexpr)
    simLRTS <- scanLRTS(X=Xsim, r=r,
                         method=method, alternative=alternative,
                         baseline=baseline,
                         case=case,
                         ...)
    sim[i] <- max(simLRTS)
  }
  pval <- mean(c(sim,obs) >= obs, na.rm=TRUE)
  names(obs) <- "maxLRTS"
  nm.alternative <- switch(alternative,
                           greater="Excess of",
                           less="Deficit of",
                           two.sided="Two-sided: excess or deficit of",
                           stop("Unknown alternative"))
  nm.alternative <- paste(nm.alternative, counted)
  result <- list(statistic = obs,
                 p.value = pval,
                 alternative = nm.alternative, 
                 method = methodname,
                 data.name = dataname)
  class(result) <- c("scan.test", "htest")
  attr(result, "obsLRTS") <- obsLRTS
  attr(result, "X") <- X
  attr(result, "r") <- r
  return(result)
}

plot.scan.test <- function(x, ..., what=c("statistic", "radius"),
                           do.window=TRUE) {
  xname <- short.deparse(substitute(x))
  what <- match.arg(what)
  Z <- as.im(x, what=what)
  do.call(plot, resolve.defaults(list(x=Z), list(...), list(main=xname)))
  if(do.window) {
    X <- attr(x, "X")
    plot(as.owin(X), add=TRUE, invert=TRUE)
  }
  invisible(NULL)
}

as.im.scan.test <- function(X, ..., what=c("statistic", "radius")) {
  Y <- attr(X, "obsLRTS")
  what <- match.arg(what)
  if(what == "radius") {
    iopt <- attr(Y, "iopt")
    r <- attr(X, "r")
    Y <- eval.im(r[iopt])
  }
  return(as.im(Y, ...))
}
