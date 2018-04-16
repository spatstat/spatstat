#
#  mincontrast.R
#
#  Functions for estimation by minimum contrast
#

##################  base ################################

mincontrast <- local({

  ## objective function (in a format that is re-usable by other code)
  contrast.objective <- function(par, objargs, ...) {
    with(objargs, {
      theo <- theoretical(par=par, rvals, ...)
      if(!is.vector(theo) || !is.numeric(theo))
        stop("theoretical function did not return a numeric vector")
      if(length(theo) != nrvals)
        stop("theoretical function did not return the correct number of values")
      if(!is.null(adjustment)) {
        theo <- adjustment$fun(theo=theo, par=par, auxdata=adjustment$auxdata)
        if(!is.vector(theo) || !is.numeric(theo))
	  stop("adjustment did not return a numeric vector")
        if(length(theo) != nrvals)
          stop("adjustment did not return the correct number of values")
      }	
      discrep <- (abs(theo^qq - obsq))^pp
      value <- mean(discrep)
      value <- min(value, .Machine$double.xmax)
      return(value)
    })
  }

  mincontrast <- function(observed, theoretical, startpar,
                          ...,
                          ctrl=list(q = 1/4, p = 2, rmin=NULL, rmax=NULL),
                          fvlab=list(label=NULL, desc="minimum contrast fit"),
                          explain=list(dataname=NULL,
                            modelname=NULL, fname=NULL),
			  adjustment=NULL) {
    verifyclass(observed, "fv")

    stopifnot(is.function(theoretical))
    if(!any("par" %in% names(formals(theoretical))))
      stop(paste("Theoretical function does not include an argument called",
                 sQuote("par")))

    ## enforce defaults
    ctrl <- resolve.defaults(ctrl, list(q = 1/4, p = 2, rmin=NULL, rmax=NULL))
    fvlab <- resolve.defaults(fvlab,
                              list(label=NULL, desc="minimum contrast fit"))
    explain <- resolve.defaults(explain,
                                list(dataname=NULL, modelname=NULL, fname=NULL))
  
    ## extract vector of r values
    argu <- fvnames(observed, ".x")
    rvals <- observed[[argu]]
    
    ## determine range of r values
    rmin <- ctrl$rmin
    rmax <- ctrl$rmax
    if(!is.null(rmin) && !is.null(rmax)) 
      stopifnot(rmin < rmax && rmin >= 0)
    else {
      alim <- attr(observed, "alim") %orifnull% range(rvals)
      if(is.null(rmax)) rmax <- alim[2]
      if(is.null(rmin)) {
        rmin <- alim[1]
        if(rmin == 0 && identical(explain$fname,"g"))
          rmin <- rmax/1e3 # avoid artefacts at zero in pcf
      }
    }
    ## extract vector of observed values of statistic
    valu <- fvnames(observed, ".y")
    obs <- observed[[valu]]
    ## restrict to [rmin, rmax]
    if(max(rvals) < rmax)
      stop(paste("rmax=", signif(rmax,4),
                 "exceeds the range of available data",
                 "= [", signif(min(rvals),4), ",", signif(max(rvals),4), "]"))
    sub <- (rvals >= rmin) & (rvals <= rmax)
    rvals <- rvals[sub]
    obs <- obs[sub]
    ## sanity clause
    if(!all(ok <- is.finite(obs))) {
      whinge <- paste("Some values of the empirical function",
                      sQuote(explain$fname),
                      "were infinite or NA.")
      iMAX <- max(which(ok))
      iMIN <- min(which(!ok)) + 1
      if(iMAX > iMIN && all(ok[iMIN:iMAX])) {
        rmin <- rvals[iMIN]
        rmax <- rvals[iMAX]
        obs   <- obs[iMIN:iMAX]
        rvals <- rvals[iMIN:iMAX]
        sub[sub] <- ok
        warning(paste(whinge,
                      "Range of r values was reset to",
                      prange(c(rmin, rmax))),
                call.=FALSE)
      } else stop(paste(whinge, "Please choose a narrower range [rmin, rmax]"),
                  call.=FALSE)
    }
    ## pack data into a list
    objargs <- list(theoretical = theoretical,
                    rvals       = rvals,
                    nrvals      = length(rvals),
                    obsq        = obs^(ctrl$q),   ## for efficiency
                    qq          = ctrl$q,
                    pp          = ctrl$p,
                    rmin        = rmin,
                    rmax        = rmax,
		    adjustment  = adjustment)
    ## go
    minimum <- optim(startpar, fn=contrast.objective, objargs=objargs, ...)
    ## if convergence failed, issue a warning 
    signalStatus(optimStatus(minimum), errors.only=TRUE)
    ## evaluate the fitted theoretical curve
    fittheo <- theoretical(minimum$par, rvals, ...)
    ## pack it up as an `fv' object
    label <- fvlab$label %orifnull% "%s[fit](r)"
    desc  <- fvlab$desc
    fitfv <- bind.fv(observed[sub, ],
                     data.frame(fit=fittheo),
                     label, desc)
    if(!is.null(adjustment)) {
      adjtheo <- adjustment$fun(theo=fittheo,
      	                        par=minimum$par,
				auxdata=adjustment$auxdata)
      fitfv <- bind.fv(fitfv,
                       data.frame(adjfit=adjtheo),
		       "%s[adjfit](r)",
		       paste("adjusted", desc))
    }				
    result <- list(par      = minimum$par,
                   fit      = fitfv,
                   opt      = minimum,
                   ctrl     = list(p=ctrl$p,q=ctrl$q,rmin=rmin,rmax=rmax),
                   info     = explain,
                   startpar = startpar,
                   objfun   = contrast.objective,
                   objargs  = objargs,
                   dotargs  = list(...))
    class(result) <- c("minconfit", class(result))
    return(result)
  }

  mincontrast
})

print.minconfit <- function(x, ...) {
  terselevel <- spatstat.options('terse')
  digits <- getOption('digits')
  ## explanatory
  cat(paste("Minimum contrast fit ",
            "(",
            "object of class ",
            dQuote("minconfit"),
            ")",
            "\n", sep=""))
  mo <- x$info$modelname
  fu <- x$info$fname
  da <- x$info$dataname
  cm <- x$covmodel
  if(!is.null(mo))
    cat("Model:", mo, fill=TRUE)
  if(!is.null(cm)) {
    ## Covariance/kernel model and nuisance parameters 
    cat("\t", cm$type, "model:", cm$model, fill=TRUE)
    margs <- cm$margs
    if(!is.null(margs)) {
      nama <- names(margs)
      tags <- ifelse(nzchar(nama), paste(nama, "="), "")
      tagvalue <- paste(tags, margs)
      splat("\t", cm$type, "parameters:",
            paste(tagvalue, collapse=", "))
    }
  }
  if(!is.null(fu) && !is.null(da))
    splat("Fitted by matching theoretical", fu, "function to", da)
  else {
    if(!is.null(fu))
      splat(" based on", fu)
    if(!is.null(da))
      splat(" fitted to", da)
  }

  if(waxlyrical('space', terselevel))
      cat("\n")
  ## Values
  splat("Internal parameters fitted by minimum contrast ($par):")
  print(x$par, ...)
  if(waxlyrical('space', terselevel))
      cat("\n")
  
  ## Handling new parameters
  isPCP <- x$isPCP %orifnull% x$internal$model!="lgcp"
  cpar <- x$clustpar
  if (!is.null(cpar)) {
    splat("Fitted",
          if(isPCP) "cluster" else "covariance",
          "parameters:")
    print(cpar, digits=digits)
  } else{
    ## Old modelpar field if necessary
    mp <- x$modelpar
    if(!is.null(mp)) {
      splat("Derived parameters of",
            if(!is.null(mo)) mo else "model",
            "($modelpar):")
      print(mp)
    }
  }
  if(length(mu <- x$mu)) {
    if(isPCP) {
      splat("Mean cluster size: ",
            if(is.numeric(mu)) paste(signif(mu, digits), "points") else
            if(is.im(mu)) "[pixel image]" else "[unknown]")
    } else {
      splat("Fitted mean of log of random intensity: ",
            if(is.numeric(mu)) signif(mu, digits) else
            if(is.im(mu)) "[pixel image]" else "[unknown]")
    }
  }
  if(waxlyrical('space', terselevel))
      cat("\n")
  ## Diagnostics
  printStatus(optimStatus(x$opt))
  ## Starting values
  if(waxlyrical('gory', terselevel)){
      cat("\n")
      splat("Starting values of parameters:")
      print(x$startpar)
      ## Algorithm parameters
      ct <- x$ctrl
      splat("Domain of integration:",
            "[",
            signif(ct$rmin,4),
            ",",
            signif(ct$rmax,4),
            "]")
      splat("Exponents:",
            "p=", paste(signif(ct$p, 3), ",",  sep=""),
            "q=", signif(ct$q,3))
  }
  invisible(NULL)
}
              

plot.minconfit <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  do.call(plot.fv,
          resolve.defaults(list(x$fit),
                           list(...),
                           list(main=xname)))
}

unitname.minconfit <- function(x) {
  unitname(x$fit)
}

"unitname<-.minconfit" <- function(x, value) {
  unitname(x$fit) <- value
  return(x)
}

as.fv.minconfit <- function(x) x$fit

######  convergence status of 'optim' object

optimStatus <- function(x, call=NULL) {
  cgce <- x$convergence
  neval <- x$counts[["function"]]
  switch(paste(cgce),
         "0" = {
           simpleMessage(
                         paste("Converged successfully after", 
                               neval, "function evaluations"),
                         call)
         },
         "1" = simpleWarning(
           paste("Iteration limit maxit was reached after",
                 neval, "function evaluations"),
           call),
         "10" = simpleWarning("Nelder-Mead simplex was degenerate", call),
         "51"= {
           simpleWarning(
                         paste("Warning message from L-BGFS-B method:",
                               sQuote(x$message)),
                         call)
         },
         "52"={
           simpleError(
                         paste("Error message from L-BGFS-B method:",
                               sQuote(x$message)),
                         call)
         },
         simpleWarning(paste("Unrecognised error code", cgce), call)
         )
}

signalStatus <- function(x, errors.only=FALSE) {
  stopifnot(inherits(x, "condition"))
  if(inherits(x, "error")) stop(x)
  if(inherits(x, "warning")) warning(x) 
  if(inherits(x, "message") && !errors.only) message(x)
  return(invisible(NULL))
}

printStatus <- function(x, errors.only=FALSE) {
  prefix <-
    if(inherits(x, "error")) "error: " else 
    if(inherits(x, "warning")) "warning: " else NULL
  if(!is.null(prefix) || !errors.only)
    cat(paste(prefix, conditionMessage(x), "\n", sep=""))
  return(invisible(NULL))
}

accumulateStatus <- function(x, stats=NULL) {
  if(is.null(stats))
    stats <- list(values=list(), frequencies=integer(0))
  if(!inherits(x, c("error", "warning", "message")))
    return(stats)
  with(stats,
       {
         same <- unlist(lapply(values, identical, y=x))
         if(any(same)) {
           i <- min(which(same))
           frequencies[i] <- frequencies[i] + 1
         } else {
           values <- append(values, list(x))
           frequencies <- c(frequencies, 1)
         }
       })
  stats <- list(values=values, frequencies=frequencies)
  return(stats)
}

printStatusList <- function(stats) {
  with(stats,
       {
         for(i in seq_along(values)) {
           printStatus(values[i])
           cat(paste("\t", paren(paste(frequencies[i], "times")), "\n"))
         }
       }
       )
  invisible(NULL)
}

  
############### applications (specific models) ##################


getdataname <- function(defaultvalue, ..., dataname=NULL) {
  if(!is.null(dataname)) dataname else defaultvalue
}
  
thomas.estK <- function(X, startpar=c(kappa=1,scale=1),
                        lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...) {

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)

  if(inherits(X, "fv")) {
    K <- X
    if(!identical(attr(K, "fname")[1], "K"))
      warning("Argument X does not appear to be a K-function")
  } else if(inherits(X, "ppp")) {
    K <- Kest(X)
    dataname <- paste("Kest(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  info <- spatstatClusterModelInfo("Thomas")
  startpar <- info$checkpar(startpar)
  theoret <- info$K
  
  result <- mincontrast(K, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)",
                          desc="minimum contrast fit of Thomas process"),
                        explain=list(dataname=dataname,
                          fname=attr(K, "fname"),
                          modelname="Thomas process"), ...)
  ## imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "sigma2")
  result$par <- par
  ## infer meaningful model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="Thomas")
  ## add new parametrisation to object
  result$clustpar <- info$checkpar(par, old=FALSE)
  return(result)
}

lgcp.estK <- function(X, startpar=c(var=1,scale=1),
                      covmodel=list(model="exponential"), 
                      lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...) {

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)
  
  if(inherits(X, "fv")) {
    K <- X
    if(!identical(attr(K, "fname")[1], "K"))
      warning("Argument X does not appear to be a K-function")
  } else if(inherits(X, "ppp")) {
    K <- Kest(X)
    dataname <- paste("Kest(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  info <- spatstatClusterModelInfo("LGCP")
  startpar <- info$checkpar(startpar)

  ## digest parameters of Covariance model and test validity
  ph <- info$parhandler
  cmodel <- do.call(ph, covmodel)
  
  theoret <- info$K

  result <- mincontrast(K, theoret, startpar,
                        ctrl=list(q=q, p=p, rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)",
                          desc="minimum contrast fit of LGCP"),
                        explain=list(dataname=dataname,
                          fname=attr(K, "fname"),
                          modelname="log-Gaussian Cox process"),
                        ...,
                        model=cmodel$model,
                        margs=cmodel$margs)
  ## imbue with meaning
  par <- result$par
  names(par) <- c("sigma2", "alpha")
  result$par <- par
  result$covmodel <- cmodel
  ## infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="lgcp")
  ## add new parametrisation to object
  result$clustpar <- info$checkpar(par, old=FALSE)
  result$clustargs <- info$checkclustargs(cmodel$margs, old=FALSE)
  return(result)
}

matclust.estK <- function(X, startpar=c(kappa=1,scale=1),
                          lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...) {

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)

  if(inherits(X, "fv")) {
    K <- X
    if(!identical(attr(K, "fname")[1], "K"))
      warning("Argument X does not appear to be a K-function")
  } else if(inherits(X, "ppp")) {
    K <- Kest(X)
    dataname <- paste("Kest(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  info <- spatstatClusterModelInfo("MatClust")
  startpar <- info$checkpar(startpar)
  theoret <- info$K
  funaux <-  info$funaux
  
  result <- mincontrast(K, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)",
                          desc="minimum contrast fit of Matern Cluster process"),
                        explain=list(dataname=dataname,
                          fname=attr(K, "fname"),
                          modelname="Matern Cluster process"),
                        ...,
                        funaux=funaux)
  ## imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "R")
  result$par <- par
  ## infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="MatClust")
  ## add new parametrisation to object
  result$clustpar <- info$checkpar(par, old=FALSE)
  return(result)
}

## versions using pcf (suggested by Jan Wild)

thomas.estpcf <- function(X, startpar=c(kappa=1,scale=1),
                          lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...,
                          pcfargs=list()){

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)

  if(inherits(X, "fv")) {
    g <- X
    if(!identical(attr(g, "fname")[1], "g"))
      warning("Argument X does not appear to be a pair correlation function")
  } else if(inherits(X, "ppp")) {
    g <- do.call(pcf.ppp, append(list(X), pcfargs))
    dataname <- paste("pcf(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  info <- spatstatClusterModelInfo("Thomas")
  startpar <- info$checkpar(startpar)
  theoret <- info$pcf
  
  ## avoid using g(0) as it may be infinite
  argu <- fvnames(g, ".x")
  rvals <- g[[argu]]
  if(rvals[1] == 0 && (is.null(rmin) || rmin == 0)) {
    rmin <- rvals[2]
  }
  result <- mincontrast(g, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(
                          label="%s[fit](r)",
                          desc="minimum contrast fit of Thomas process"),
                        explain=list(
                          dataname=dataname,
                          fname=attr(g, "fname"),
                          modelname="Thomas process"), ...)
  ## imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "sigma2")
  result$par <- par
  ## infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="Thomas")
  ## add new parametrisation to object
  result$clustpar <- info$checkpar(par, old=FALSE)
  return(result)
}

matclust.estpcf <- function(X, startpar=c(kappa=1,scale=1),
                            lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...,
                            pcfargs=list()){

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)

  if(inherits(X, "fv")) {
    g <- X
    if(!identical(attr(g, "fname")[1], "g"))
      warning("Argument X does not appear to be a pair correlation function")
  } else if(inherits(X, "ppp")) {
    g <- do.call(pcf.ppp, append(list(X), pcfargs))
    dataname <- paste("pcf(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  info <- spatstatClusterModelInfo("MatClust")
  startpar <- info$checkpar(startpar)
  theoret <- info$pcf
  funaux <-  info$funaux
  
  ## avoid using g(0) as it may be infinite
  argu <- fvnames(g, ".x")
  rvals <- g[[argu]]
  if(rvals[1] == 0 && (is.null(rmin) || rmin == 0)) {
    rmin <- rvals[2]
  }
  result <- mincontrast(g, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)",
                          desc="minimum contrast fit of Matern Cluster process"),
                        explain=list(dataname=dataname,
                          fname=attr(g, "fname"),
                          modelname="Matern Cluster process"),
                        ...,
                        funaux=funaux)
  ## imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "R")
  result$par <- par
  ## infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="MatClust")
  ## add new parametrisation to object
  result$clustpar <- info$checkpar(par, old=FALSE)
  return(result)
}

lgcp.estpcf <- function(X, startpar=c(var=1,scale=1),
                      covmodel=list(model="exponential"), 
                        lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...,
                        pcfargs=list()) {
  
  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)
  
  if(inherits(X, "fv")) {
    g <- X
    if(!identical(attr(g, "fname")[1], "g"))
      warning("Argument X does not appear to be a pair correlation function")
  } else if(inherits(X, "ppp")) {
    g <- do.call(pcf.ppp, append(list(X), pcfargs))
    dataname <- paste("pcf(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  info <- spatstatClusterModelInfo("LGCP")
  startpar <- info$checkpar(startpar)

  ## digest parameters of Covariance model and test validity
  ph <- info$parhandler
  cmodel <- do.call(ph, covmodel)
  
  theoret <- info$pcf
  
  result <- mincontrast(g, theoret, startpar,
                        ctrl=list(q=q, p=p, rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)",
                          desc="minimum contrast fit of LGCP"),
                        explain=list(dataname=dataname,
                          fname=attr(g, "fname"),
                          modelname="log-Gaussian Cox process"),
                        ...,
                        model=cmodel$model,
                        margs=cmodel$margs)
  ## imbue with meaning
  par <- result$par
  names(par) <- c("sigma2", "alpha")
  result$par <- par
  result$covmodel <- cmodel
  ## infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="lgcp")
  ## add new parametrisation to object
  result$clustpar <- info$checkpar(par, old=FALSE)
  result$clustargs <- info$checkclustargs(cmodel$margs, old=FALSE)
  return(result)
}


cauchy.estK <- function(X, startpar=c(kappa=1,scale=1),
                        lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...) {

## omega: scale parameter of Cauchy kernel function
## eta: scale parameter of Cauchy pair correlation function
## eta = 2 * omega

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)

  if(inherits(X, "fv")) {
    K <- X
    if(!identical(attr(K, "fname")[1], "K"))
      warning("Argument X does not appear to be a K-function")
  } else if(inherits(X, "ppp")) {
    K <- Kest(X)
    dataname <- paste("Kest(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  info <- spatstatClusterModelInfo("Cauchy")
  startpar <- info$checkpar(startpar)
  theoret <- info$K

  desc <- "minimum contrast fit of Neyman-Scott process with Cauchy kernel"
  result <- mincontrast(K, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)", desc=desc),
                        explain=list(dataname=dataname,
                          fname=attr(K, "fname"),
                          modelname="Cauchy process"), ...)
  ## imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "eta2")
  result$par <- par
  ## infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="Cauchy")
  ## add new parametrisation to object
  result$clustpar <- info$checkpar(par, old=FALSE)
  return(result)
}


cauchy.estpcf <- function(X, startpar=c(kappa=1,scale=1),
                          lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...,
                          pcfargs=list()) {

## omega: scale parameter of Cauchy kernel function
## eta: scale parameter of Cauchy pair correlation function
## eta = 2 * omega

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)

  if(inherits(X, "fv")) {
    g <- X
    if(!identical(attr(g, "fname")[1], "g"))
      warning("Argument X does not appear to be a pair correlation function")
  } else if(inherits(X, "ppp")) {
    g <- do.call(pcf.ppp, append(list(X), pcfargs))
    dataname <- paste("pcf(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  info <- spatstatClusterModelInfo("Cauchy")
  startpar <- info$checkpar(startpar)
  theoret <- info$pcf

  ## avoid using g(0) as it may be infinite
  argu <- fvnames(g, ".x")
  rvals <- g[[argu]]
  if(rvals[1] == 0 && (is.null(rmin) || rmin == 0)) {
    rmin <- rvals[2]
  }
  
  desc <- "minimum contrast fit of Neyman-Scott process with Cauchy kernel"
  result <- mincontrast(g, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)", desc=desc),
                        explain=list(dataname=dataname,
                          fname=attr(g, "fname"),
                          modelname="Cauchy process"), ...)
  ## imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "eta2")
  result$par <- par
  ## infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="Cauchy")
  ## add new parametrisation to object
  result$clustpar <- info$checkpar(par, old=FALSE)
  return(result)
}

## user-callable
resolve.vargamma.shape <- function(..., nu.ker=NULL, nu.pcf=NULL, default = FALSE) {
  if(is.null(nu.ker) && is.null(nu.pcf)){
    if(!default)
        stop("Must specify either nu.ker or nu.pcf", call.=FALSE)
    nu.ker <- -1/4
  }
  if(!is.null(nu.ker) && !is.null(nu.pcf))
    stop("Only one of nu.ker and nu.pcf should be specified",
         call.=FALSE)
  if(!is.null(nu.ker)) {
    check.1.real(nu.ker)
    stopifnot(nu.ker > -1/2)
    nu.pcf <- 2 * nu.ker + 1
  } else {
    check.1.real(nu.pcf)
    stopifnot(nu.pcf > 0)
    nu.ker <- (nu.pcf - 1)/2
  }
  return(list(nu.ker=nu.ker, nu.pcf=nu.pcf))
}

vargamma.estK <- function(X, startpar=c(kappa=1,scale=1), nu = -1/4,
                          lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL,
                          ...) {

## nu.ker: smoothness parameter of Variance Gamma kernel function
## omega: scale parameter of kernel function
## nu.pcf: smoothness parameter of Variance Gamma pair correlation function
## eta: scale parameter of Variance Gamma pair correlation function
## nu.pcf = 2 * nu.ker + 1    and    eta = omega

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)
  
  if(inherits(X, "fv")) {
    K <- X
    if(!identical(attr(K, "fname")[1], "K"))
      warning("Argument X does not appear to be a K-function")
  } else if(inherits(X, "ppp")) {
    K <- Kest(X)
    dataname <- paste("Kest(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  ## Catch old nu.ker/nu.pcf syntax and resolve nu-value.
  dots <- list(...)
  if(missing(nu)){
      nu <- resolve.vargamma.shape(nu.ker=dots$nu.ker, nu.pcf=dots$nu.pcf, default = TRUE)$nu.ker
  }
  check.1.real(nu)
  stopifnot(nu > -1/2)

  info <- spatstatClusterModelInfo("VarGamma")
  startpar <- info$checkpar(startpar)
  theoret <- info$K
  
  ## test validity of parameter nu and digest
  ph <- info$parhandler
  cmodel <- ph(nu.ker=nu)
  margs <- cmodel$margs

  desc <- "minimum contrast fit of Neyman-Scott process with Variance Gamma kernel"
  result <- mincontrast(K, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)", desc=desc),
                        explain=list(dataname=dataname,
                          fname=attr(K, "fname"),
                          modelname="Variance Gamma process"),
                        margs=margs, ...)
  ## imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "eta")
  result$par <- par
  result$covmodel <- cmodel
  ## infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="VarGamma")
  ## add new parametrisation to object
  result$clustpar <- info$checkpar(par, old=FALSE)
  result$clustargs <- info$checkclustargs(cmodel$margs, old=FALSE)
  return(result)
}


vargamma.estpcf <- function(X, startpar=c(kappa=1,scale=1), nu=-1/4, 
                            lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, 
                            ..., pcfargs=list()) {

## nu.ker: smoothness parameter of Variance Gamma kernel function
## omega: scale parameter of kernel function
## nu.pcf: smoothness parameter of Variance Gamma pair correlation function
## eta: scale parameter of Variance Gamma pair correlation function
## nu.pcf = 2 * nu.ker + 1    and    eta = omega

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)

  if(inherits(X, "fv")) {
    g <- X
    if(!identical(attr(g, "fname")[1], "g"))
      warning("Argument X does not appear to be a pair correlation function")
  } else if(inherits(X, "ppp")) {
    g <- do.call(pcf.ppp, append(list(X), pcfargs))
    dataname <- paste("pcf(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
      stop("Unrecognised format for argument X")
  
  ## Catch old nu.ker/nu.pcf syntax and resolve nu-value.
  dots <- list(...)
  if(missing(nu)){
      ## nutmp <- try(resolve.vargamma.shape(nu.ker=dots$nu.ker, nu.pcf=dots$nu.pcf)$nu.ker, silent=TRUE)
      ## if(!inherits(nutmp, "try-error")) nu <- nutmp
      nu <- resolve.vargamma.shape(nu.ker=dots$nu.ker, nu.pcf=dots$nu.pcf, default = TRUE)$nu.ker
  }
  check.1.real(nu)
  stopifnot(nu > -1/2)

  info <- spatstatClusterModelInfo("VarGamma")
  startpar <- info$checkpar(startpar)
  theoret <- info$pcf

  ## test validity of parameter nu and digest 
  ph <- info$parhandler
  cmodel <- ph(nu.ker=nu)
  margs <- cmodel$margs
  
  ## avoid using g(0) as it may be infinite
  argu <- fvnames(g, ".x")
  rvals <- g[[argu]]
  if(rvals[1] == 0 && (is.null(rmin) || rmin == 0)) {
    rmin <- rvals[2]
  }
  
  desc <- "minimum contrast fit of Neyman-Scott process with Variance Gamma kernel"
  result <- mincontrast(g, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)", desc=desc),
                        explain=list(dataname=dataname,
                          fname=attr(g, "fname"),
                          modelname="Variance Gamma process"),
                        margs=margs,
                        ...)
  ## imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "eta")
  result$par <- par
  result$covmodel <- cmodel
  ## infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="VarGamma")
  ## add new parametrisation to object
  result$clustpar <- info$checkpar(par, old=FALSE)
  result$clustargs <- info$checkclustargs(cmodel$margs, old=FALSE)
  return(result)
}

