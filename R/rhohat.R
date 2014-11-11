#
#  rhohat.R
#
#  $Revision: 1.61 $  $Date: 2014/11/10 11:18:05 $
#
#  Non-parametric estimation of a transformation rho(z) determining
#  the intensity function lambda(u) of a point process in terms of a
#  spatial covariate Z(u) through lambda(u) = rho(Z(u)).
#  More generally allows offsets etc.

rhohat <- function(object, covariate, ...) {
  UseMethod("rhohat")
}

rhohat.ppp <- rhohat.quad <- 
  function(object, covariate, ...,
           baseline=NULL, weights=NULL,
           method=c("ratio", "reweight", "transform"),
           horvitz=FALSE,
           smoother=c("kernel", "local"),
           dimyx=NULL, eps=NULL,
           n=512, bw="nrd0", adjust=1, from=NULL, to=NULL, 
           bwref=bw, covname, confidence=0.95) {
  callstring <- short.deparse(sys.call())
  smoother <- match.arg(smoother)
  method <- match.arg(method)
  if(missing(covname)) 
    covname <- sensiblevarname(short.deparse(substitute(covariate)), "X")
  if(is.null(adjust))
    adjust <- 1
  # validate model
  if(is.null(baseline)) {
    model <- ppm(object ~1)
    reference <- "Lebesgue"
  } else {
    model <- ppm(object ~ offset(log(baseline)))
    reference <- "baseline"
  } 
  modelcall <- NULL

  if(is.character(covariate) && length(covariate) == 1) {
    covname <- covariate
    switch(covname,
           x={
             covariate <- function(x,y) { x }
           }, 
           y={
             covariate <- function(x,y) { y }
           },
           stop("Unrecognised covariate name")
         )
    covunits <- unitname(data.ppm(model))
  } else {
    covunits <- NULL
  }

  areaW <- area(Window(data.ppm(model)))
  
  rhohatEngine(model, covariate, reference, areaW, ...,
               weights=weights,
               method=method,
               horvitz=horvitz,
               smoother=smoother,
               resolution=list(dimyx=dimyx, eps=eps),
               n=n, bw=bw, adjust=adjust, from=from, to=to,
               bwref=bwref, covname=covname, covunits=covunits,
               confidence=confidence,
               modelcall=modelcall, callstring=callstring)
}

rhohat.ppm <- function(object, covariate, ...,
                       weights=NULL,
                       method=c("ratio", "reweight", "transform"),
                       horvitz=FALSE,
                       smoother=c("kernel", "local"),
                       dimyx=NULL, eps=NULL,
                       n=512, bw="nrd0", adjust=1, from=NULL, to=NULL, 
                       bwref=bw, covname, confidence=0.95) {
  callstring <- short.deparse(sys.call())
  smoother <- match.arg(smoother)
  method <- match.arg(method)
  if(missing(covname)) 
    covname <- sensiblevarname(short.deparse(substitute(covariate)), "X")
  if(is.null(adjust))
    adjust <- 1

  if("baseline" %in% names(list(...)))
    warning("Argument 'baseline' ignored: not available for rhohat.ppm")

  ## validate model
  model <- object
  reference <- "model"
  modelcall <- model$call

  if(is.character(covariate) && length(covariate) == 1) {
    covname <- covariate
    switch(covname,
           x={
             covariate <- function(x,y) { x }
           }, 
           y={
             covariate <- function(x,y) { y }
           },
           stop("Unrecognised covariate name")
         )
    covunits <- unitname(data.ppm(model))
  } else {
    covunits <- NULL
  }

  areaW <- area(Window(data.ppm(model)))
  
  rhohatEngine(model, covariate, reference, areaW, ...,
               weights=weights,
               method=method,
               horvitz=horvitz,
               smoother=smoother,
               resolution=list(dimyx=dimyx, eps=eps),
               n=n, bw=bw, adjust=adjust, from=from, to=to,
               bwref=bwref, covname=covname, covunits=covunits,
               confidence=confidence,
               modelcall=modelcall, callstring=callstring)
}

rhohat.lpp <- rhohat.lppm <- 
  function(object, covariate, ...,
           weights=NULL,
           method=c("ratio", "reweight", "transform"),
           horvitz=FALSE,
           smoother=c("kernel", "local"),
           nd=1000, eps=NULL,
           n=512, bw="nrd0", adjust=1, from=NULL, to=NULL, 
           bwref=bw, covname, confidence=0.95) {
  callstring <- short.deparse(sys.call())
  smoother <- match.arg(smoother)
  method <- match.arg(method)
  if(missing(covname)) 
    covname <- sensiblevarname(short.deparse(substitute(covariate)), "X")
  if(is.null(adjust))
    adjust <- 1
  # validate model
  if(is.lpp(object)) {
    X <- object
    model <- lppm(object, ~1)
    reference <- "Lebesgue"
    modelcall <- NULL
  } else if(inherits(object, "lppm")) {
    model <- object
    X <- model$X
    reference <- "model"
    modelcall <- model$call
  } else stop("object should be of class lpp or lppm")
  
  if("baseline" %in% names(list(...)))
    warning("Argument 'baseline' ignored: not available for ",
            if(is.lpp(object)) "rhohat.lpp" else "rhohat.lppm")

  if(is.character(covariate) && length(covariate) == 1) {
    covname <- covariate
    switch(covname,
           x={
             covariate <- function(x,y) { x }
           }, 
           y={
             covariate <- function(x,y) { y }
           },
           stop("Unrecognised covariate name")
         )
    covunits <- unitname(X)
  } else {
    covunits <- NULL
  }

  totlen <- sum(lengths.psp(as.psp(as.linnet(X))))
  
  rhohatEngine(model, covariate, reference, totlen, ...,
               weights=weights,
               method=method,
               horvitz=horvitz,
               smoother=smoother,
               resolution=list(nd=nd, eps=eps),
               n=n, bw=bw, adjust=adjust, from=from, to=to,
               bwref=bwref, covname=covname, covunits=covunits,
               confidence=confidence,
               modelcall=modelcall, callstring=callstring)
}

rhohatEngine <- function(model, covariate,
                         reference=c("Lebesgue", "model", "baseline"),
                         volume,
                         ...,
                         weights=NULL,
                         method=c("ratio", "reweight", "transform"),
                         horvitz=FALSE,
                         smoother=c("kernel", "local"),
                         resolution=list(), 
                         n=512, bw="nrd0", adjust=1, from=NULL, to=NULL, 
                         bwref=bw, covname, covunits=NULL, confidence=0.95,
                         modelcall=NULL, callstring="rhohat") {
  reference <- match.arg(reference)
  # evaluate the covariate at data points and at pixels
  stuff <- do.call("evalCovar",
                   append(list(model, covariate), resolution))
  # unpack
#  info   <- stuff$info
  values <- stuff$values
  # values at each data point
  ZX      <- values$ZX
  lambdaX <- if(horvitz) fitted(model, dataonly=TRUE) else NULL
  # values at each pixel
  Zimage  <- values$Zimage
  Zvalues <- values$Zvalues
  lambda  <- values$lambda
  ## weights
  if(!is.null(weights)) {
    X <- data.ppm(model)
    if(is.im(weights)) 
      weights <- safelookup(weights, X)
    else if(is.function(weights))
      weights <- weights(X$x, X$y)
    else if(is.numeric(weights) && is.vector(as.numeric(weights))) 
      check.nvector(weights, npoints(X))
    else stop(paste(sQuote("weights"),
                    "should be a vector, a pixel image, or a function"))
  }
  # normalising constants
  denom <- volume * (if(reference == "Lebesgue" || horvitz) 1 else mean(lambda))
  # info 
  savestuff <- list(reference  = reference,
                    horvitz    = horvitz,
                    Zimage     = Zimage)
  # calculate rho-hat
  result <- rhohatCalc(ZX, Zvalues, lambda, denom,
                       ...,
                       weights=weights,
                       lambdaX=lambdaX,
                       method=method,
                       horvitz=horvitz,
                       smoother=smoother,
                       n=n, bw=bw, adjust=adjust, from=from, to=to,
                       bwref=bwref, covname=covname, confidence=confidence,
                       covunits=covunits,
                       modelcall=modelcall, callstring=callstring,
                       savestuff=savestuff)
  return(result)
}


# basic calculation of rhohat from covariate values

rhohatCalc <- local({
  
  interpolate <- function(x,y) {
    if(inherits(x, "density") && missing(y))
      approxfun(x$x, x$y, rule=2)
    else 
      approxfun(x, y, rule=2)
  }

  ## note: this function normalises the weights, like density.default
  LocfitRaw <- function(x, ..., weights=NULL) {
    if(is.null(weights)) weights <- 1
    do.call.matched("locfit.raw", append(list(x=x, weights=weights), list(...)))
  }

  varlog <- function(obj,xx) {
    ## variance of log f-hat
    stopifnot(inherits(obj, "locfit"))
    if(!identical(obj$trans, exp))
      stop("internal error: locfit object does not have log link")
    ## the following call should have band="local" but that produces NaN's
    pred <- predict(obj, newdata=xx,
                    se.fit=TRUE, what="coef")
    se <- pred$se.fit
    return(se^2)
  }

  rhohatCalc <- function(ZX, Zvalues, lambda, denom, ...,
                         weights=NULL, lambdaX=NULL,
                         method=c("ratio", "reweight", "transform"),
                         horvitz=FALSE, 
                         smoother=c("kernel", "local"),
                         n=512, bw="nrd0", adjust=1, from=NULL, to=NULL, 
                         bwref=bw, covname, confidence=0.95,
                         covunits = NULL, modelcall=NULL, callstring=NULL,
                         savestuff=list()) {
    method <- match.arg(method)
    smoother <- match.arg(smoother)
    ## check availability of locfit package
    if(smoother == "local" && !require(locfit, quietly=TRUE)) {
      warning(paste("In", paste(dQuote(callstring), ":", sep=""),
                    "package", sQuote("locfit"), "is not available;",
                    "unable to perform local likelihood smoothing;",
                    "using kernel smoothing instead"),
              call.=FALSE)
      smoother <- "kernel"
    }
    ## validate
    stopifnot(is.numeric(ZX))
    stopifnot(is.numeric(Zvalues))
    stopifnot(is.numeric(lambda))
    stopifnot(length(lambda) == length(Zvalues))
    stopifnot(all(is.finite(lambda))) 
    check.1.real(denom)
    ## 
    if(horvitz) {
      ## data points will be weighted by reciprocal of model intensity
      weights <- (weights %orifnull% 1)/lambdaX
    }
    ## normalising constants
    nX   <- if(is.null(weights)) length(ZX) else sum(weights)
    kappahat <- nX/denom
    ## limits
    Zrange <- range(ZX, Zvalues)
    if(is.null(from)) from <- Zrange[1] 
    if(is.null(to))   to   <- Zrange[2]
    if(from > Zrange[1] || to < Zrange[2])
      stop("Interval [from, to] = ", prange(c(from,to)), 
           "does not contain the range of data values =", prange(Zrange))
    ## critical constant for CI's
    crit <- qnorm((1+confidence)/2)
    percentage <- paste(round(100 * confidence), "%%", sep="")
    CIblurb <- paste("pointwise", percentage, "confidence interval")
    ## estimate densities   
    if(smoother == "kernel") {
      ## ............... kernel smoothing ......................
      ## reference density (normalised) for calculation
      ghat <- density(Zvalues,weights=if(horvitz) NULL else lambda/sum(lambda),
                      bw=bwref,adjust=adjust,n=n,from=from,to=to, ...)
      xxx <- ghat$x
      ghatfun <- interpolate(ghat)
      ## relative density
      switch(method,
             ratio={
               ## compute ratio of smoothed densities
               fhat <- unnormdensity(ZX,weights=weights,
                                     bw=bw,adjust=adjust,
                                     n=n,from=from, to=to, ...)
               fhatfun <- interpolate(fhat)
               Ghat.xxx <- denom * ghatfun(xxx)
               yyy <- fhatfun(xxx)/Ghat.xxx
               ## compute variance approximation
               sigma <- fhat$bw
               weights2 <- if(is.null(weights)) NULL else weights^2
               fstar <- unnormdensity(ZX,weights=weights2,
                                      bw=bw,adjust=adjust/sqrt(2),
                                      n=n,from=from, to=to, ...)
               fstarfun <- interpolate(fstar)
               const <- 1/(2 * sigma * sqrt(pi))
               vvv  <- const * fstarfun(xxx)/Ghat.xxx^2
             },
             reweight={
               ## weight Z values by reciprocal of reference
               wt <- (weights %orifnull% 1)/(denom * ghatfun(ZX))
               rhat <- unnormdensity(ZX, weights=wt, bw=bw,adjust=adjust,
                                     n=n,from=from, to=to, ...)
               rhatfun <- interpolate(rhat)
               yyy <- rhatfun(xxx)
               ## compute variance approximation
               sigma <- rhat$bw
               rongstar <- unnormdensity(ZX, weights=wt^2,
                                         bw=bw,adjust=adjust/sqrt(2),
                                         n=n,from=from, to=to, ...)
               rongstarfun <- interpolate(rongstar)
               const <- 1/(2 * sigma * sqrt(pi))
               vvv  <- const * rongstarfun(xxx)
             },
             transform={
               ## probability integral transform
               Gfun <- interpolate(ghat$x, cumsum(ghat$y)/sum(ghat$y))
               GZX <- Gfun(ZX)
               ## smooth density on [0,1]
               qhat <- unnormdensity(GZX,weights=weights,
                                     bw=bw,adjust=adjust,
                                     n=n, from=0, to=1, ...)
               qhatfun <- interpolate(qhat)
               ## edge effect correction
               one <- density(seq(from=0,to=1,length.out=512),
                              bw=qhat$bw, adjust=1,
                              n=n,from=0, to=1, ...)
               onefun <- interpolate(one)
               ## apply to transformed values
               Gxxx <- Gfun(xxx)
               Dxxx <- denom * onefun(Gxxx)
               yyy <- qhatfun(Gxxx)/Dxxx
               ## compute variance approximation
               sigma <- qhat$bw
               weights2 <- if(is.null(weights)) NULL else weights^2
               qstar <- unnormdensity(GZX,weights=weights2,
                                      bw=bw,adjust=adjust/sqrt(2),
                                      n=n,from=0, to=1, ...)
               qstarfun <- interpolate(qstar)
               const <- 1/(2 * sigma * sqrt(pi))
               vvv  <- const * qstarfun(Gxxx)/Dxxx^2
             })
      vvvname <- "Variance of estimator"
      vvvlabel <- paste("bold(Var)~hat(%s)", paren(covname), sep="")
      sd <- sqrt(vvv)
      hi <- yyy + crit * sd
      lo <- yyy - crit * sd
    } else {
      ## .................. local likelihood smoothing .......................
      xlim <- c(from, to)
      xxx <- seq(from, to, length=n)
      ## reference density
      ghat <- LocfitRaw(Zvalues,
                        weights=if(horvitz) NULL else lambda,
                        xlim=xlim, ...)
      ggg <- predict(ghat, xxx)
      ## relative density
      switch(method,
             ratio={
               ## compute ratio of smoothed densities
               fhat <- LocfitRaw(ZX, weights=weights, xlim=xlim, ...)
               fff <- predict(fhat, xxx)
               yyy <- kappahat * fff/ggg
               ## compute approximation to variance of log rho-hat
               varlogN <- 1/nX
               vvv <- varlog(fhat, xxx) + varlogN
             },
             reweight={
               ## weight Z values by reciprocal of reference
               wt <- (weights %orifnull% 1)/(denom * predict(ghat,ZX))
               sumwt <- sum(wt)
               rhat <- LocfitRaw(ZX, weights=wt, xlim=xlim, ...)
               rrr <- predict(rhat, xxx)
               yyy <- sumwt * rrr
               ## compute approximation to variance of log rho-hat
               varsumwt <- mean(yyy /(denom * ggg)) * diff(xlim)
               varlogsumwt <- varsumwt/sumwt^2
               vvv <- varlog(rhat, xxx) + varlogsumwt
             },
             transform={
               ## probability integral transform
               Gfun <- approxfun(xxx, cumsum(ggg)/sum(ggg), rule=2)
               GZX <- Gfun(ZX)
               ## smooth density on [0,1], end effect corrected
               qhat <- LocfitRaw(GZX, weights=weights, xlim=c(0,1), ...)
               ## apply to transformed values
               Gxxx <- Gfun(xxx)
               qqq <- predict(qhat, Gxxx)
               yyy <- kappahat * qqq
               ## compute approximation to variance of log rho-hat
               varlogN <- 1/nX
               vvv <- varlog(qhat, Gxxx) + varlogN
             })
      vvvname <- "Variance of log of estimator"
      vvvlabel <- paste("bold(Var)~log(hat(%s)", paren(covname), ")", sep="")
      sss <- exp(crit * sqrt(vvv))
      hi <- yyy * sss
      lo <- yyy / sss
    }
    ## pack into fv object
    df <- data.frame(xxx=xxx, rho=yyy, var=vvv, hi=hi, lo=lo)
    names(df)[1] <- covname
    desc <- c(paste("covariate", covname),
              "Estimated intensity",
              vvvname,
              paste("Upper limit of", CIblurb),
              paste("Lower limit of", CIblurb))
    rslt <- fv(df,
               argu=covname,
               ylab=substitute(rho(X), list(X=as.name(covname))),
               valu="rho",
               fmla= as.formula(paste(". ~ ", covname)),
               alim=range(ZX),
               labl=c(covname,
                 paste("hat(%s)", paren(covname), sep=""),
                 vvvlabel,
                 paste("%s[hi]", paren(covname), sep=""),
                 paste("%s[lo]", paren(covname), sep="")),
               desc=desc,
               unitname=covunits,
               fname="rho",
               yexp=substitute(rho(X), list(X=as.name(covname))))
    attr(rslt, "dotnames") <- c("rho", "hi", "lo")
    ## pack up
    class(rslt) <- c("rhohat", class(rslt))
    ## add info
    stuff <- 
      list(modelcall  = modelcall, 
           callstring = callstring,
           sigma      = switch(smoother, kernel=sigma, local=NULL),
           covname    = paste(covname, collapse=""),
           ZX         = ZX,
           lambda     = lambda,
           method     = method,
           smoother   = smoother)
    attr(rslt, "stuff") <- append(stuff, savestuff)
    return(rslt)
  }
  
  rhohatCalc
})

## ........... end of 'rhohatCalc' .................................


print.rhohat <- function(x, ...) {
  s <- attr(x, "stuff")
  splat("Intensity function estimate (class rhohat)",
        "for the covariate", s$covname)
  switch(s$reference,
         Lebesgue=splat("Function values are absolute intensities"),
         baseline=splat("Function values are relative to baseline"),
         model={
           splat("Function values are relative to fitted model")
           print(s$modelcall)
         })
  cat("Estimation method: ")
  switch(s$method,
         ratio={
           splat("ratio of fixed-bandwidth kernel smoothers")
         },
         reweight={
           splat("fixed-bandwidth kernel smoother of weighted data")
         },
         transform={
           splat("probability integral transform,",
                 "edge-corrected fixed bandwidth kernel smoothing",
                 "on [0,1]")
         },
         cat("UNKNOWN\n"))
  if(identical(s$horvitz, TRUE))
    splat("\twith Horvitz-Thompson weight")
  cat("Smoother: ")
  switch(s$smoother,
         kernel={
           splat("Kernel density estimator")
           splat("Actual smoothing bandwidth sigma = ",
                 signif(s$sigma,5))
         },
         local ={ splat("Local likelihood density estimator") }
         )
  splat("Call:", s$callstring)

  NextMethod("print")
}

plot.rhohat <- function(x, ..., do.rug=TRUE) {
  xname <- short.deparse(substitute(x))
  s <- attr(x, "stuff")
  covname <- s$covname
  asked.rug <- !missing(do.rug) && identical(rug, TRUE)
  out <- do.call("plot.fv",
                 resolve.defaults(list(x=x), list(...),
                                  list(main=xname, shade=c("hi", "lo"))))
  if(identical(list(...)$limitsonly, TRUE))
    return(out)
  if(do.rug) {
    rugx <- ZX <- s$ZX
    # check whether it's the default plot
    argh <- list(...)
    isfo <- unlist(lapply(argh, inherits, what="formula"))
    if(any(isfo)) {
      # a plot formula was given; inspect RHS
      fmla <- argh[[min(which(isfo))]]
      rhs <- rhs.of.formula(fmla)
      vars <- variablesinformula(rhs)
      vars <- vars[vars %in% c(colnames(x), ".x", ".y")]
      if(length(vars) == 1 && vars %in% c(covname, ".x")) {
        # expression in terms of covariate
        rhstr <- as.character(rhs)[2]
        dat <- list(ZX)
        names(dat) <- vars[1]
        rugx <- as.numeric(eval(parse(text=rhstr), dat))
      } else {
        if(asked.rug) warning("Unable to add rug plot")
        rugx <- NULL
      }
    } 
    if(!is.null(rugx)) {
      # restrict to x limits, if given
      if(!is.null(xlim <- list(...)$xlim))
        rugx <- rugx[rugx >= xlim[1] & rugx <= xlim[2]]
      # finally plot the rug
      if(length(rugx) > 0)
        rug(rugx)
    }
  }
  invisible(NULL)
}

predict.rhohat <- function(object, ..., relative=FALSE) {
  if(length(list(...)) > 0)
    warning("Additional arguments ignored in predict.rhohat")
  # extract info
  s <- attr(object, "stuff")
  reference <- s$reference
  # convert to (linearly interpolated) function 
  x <- with(object, .x)
  y <- with(object, .y)
  rho <- approxfun(x, y, rule=2)
  # extract image of covariate
  Z <- s$Zimage
  # apply rho to Z
  Y <- eval.im(rho(Z))
  # adjust to reference baseline
  if(reference != "Lebesgue" && !relative) {
    lambda <- s$lambda
    Y <- eval.im(Y * lambda)
  }
  return(Y)
}

as.function.rhohat <- function(x, ..., value=".y", extrapolate=TRUE) {
  NextMethod("as.function")
}
