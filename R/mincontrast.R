#
#  mincontrast.R
#
#  Functions for estimation by minimum contrast
#

##################  base ################################

mincontrast <- function(observed, theoretical, startpar,
                        ...,
                        ctrl=list(q = 1/4, p = 2, rmin=NULL, rmax=NULL),
                        fvlab=list(label=NULL, desc="minimum contrast fit"),
                        explain=list(dataname=NULL, modelname=NULL, fname=NULL)) {
  verifyclass(observed, "fv")
  stopifnot(is.function(theoretical))
  if(!any("par" %in% names(formals(theoretical))))
    stop(paste("Theoretical function does not include an argument called",
               sQuote("par")))

  # enforce defaults
  ctrl <- resolve.defaults(ctrl, list(q = 1/4, p = 2, rmin=NULL, rmax=NULL))
  fvlab <- resolve.defaults(fvlab,
                            list(label=NULL, desc="minimum contrast fit"))
  explain <- resolve.defaults(explain,
                              list(dataname=NULL, modelname=NULL, fname=NULL))
  
  # determine range of r values
  rmin <- ctrl$rmin
  rmax <- ctrl$rmax
  if(!is.null(rmin) && !is.null(rmax)) 
    stopifnot(rmin < rmax && rmin >= 0)
  else {
    alim <- attr(observed, "alim")
    if(is.null(rmin)) rmin <- alim[1]
    if(is.null(rmax)) rmax <- alim[2]
  }
  # extract vector of r values
  argu <- fvnames(observed, ".x")
  rvals <- observed[[argu]]
  # extract vector of observed values of statistic
  valu <- fvnames(observed, ".y")
  obs <- observed[[valu]]
  # restrict to [rmin, rmax]
  if(max(rvals) < rmax)
    stop(paste("rmax=", signif(rmax,4),
               "exceeds the range of available data",
               "= [", signif(min(rvals),4), ",", signif(max(rvals),4), "]"))
  sub <- (rvals >= rmin) & (rvals <= rmax)
  rvals <- rvals[sub]
  obs <- obs[sub]
  # sanity clause
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
  # for efficiency
  obsq <- obs^(ctrl$q)
  # define objective function
  objective <- function(par, obsq, theoretical, rvals, qq, pp, rmin, rmax, ...) {
    theo <- theoretical(par=par, rvals, ...)
    if(!is.vector(theo) || !is.numeric(theo))
      stop("theoretical function did not return a numeric vector")
    if(length(theo) != length(obs))
      stop("theoretical function did not return the correct number of values")
    discrep <- (abs(theo^qq - obsq))^pp
    return(sum(discrep))
  }
  # go
  minimum <- optim(startpar, fn=objective,
                   obsq=obsq, theoretical=theoretical,
                   rvals=rvals,
                   qq=ctrl$q, pp=ctrl$p, rmin=rmin, rmax=rmax, ...)
  # if convergence failed, issue a warning 
  signalStatus(optimStatus(minimum), errors.only=TRUE)
  # evaluate the fitted theoretical curve
  fittheo <- theoretical(minimum$par, rvals, ...)
  # pack it up as an `fv' object
  label <- fvlab$label
  desc  <- fvlab$desc
  if(is.null(label))
    label <- paste("fit(", argu, ")", collapse="")
  fitfv <- bind.fv(observed[sub, ],
                   data.frame(fit=fittheo),
                   label, desc)
  result <- list(par=minimum$par,
                 fit=fitfv,
                 opt=minimum,
                 ctrl=list(p=ctrl$p,q=ctrl$q,rmin=rmin,rmax=rmax),
                 info=explain,
                 startpar=startpar
                 )
  
  class(result) <- c("minconfit", class(result))
  return(result)
}

print.minconfit <- function(x, ...) {
  # explanatory
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
    cat(paste("Model:", mo, "\n"))
  if(!is.null(cm)) {
    # Covariance/kernel model and nuisance parameters 
    cat(paste("\t", cm$type, "model:", cm$model, "\n"))
    margs <- cm$margs
    if(!is.null(margs)) {
      nama <- names(margs)
      tags <- ifelse(nzchar(nama), paste(nama, "="), "")
      tagvalue <- paste(tags, margs)
      cat(paste("\t", cm$type, "parameters:",
                paste(tagvalue, collapse=", "),
                "\n"))
    }
  }
  if(!is.null(fu) && !is.null(da))
    cat(paste("Fitted by matching theoretical", fu, "function to", da))
  else {
    if(!is.null(fu))
      cat(paste(" based on", fu))
    if(!is.null(da))
      cat(paste(" fitted to", da))
  }
  cat("\n")
  # Values
  cat("Parameters fitted by minimum contrast ($par):\n")
  print(x$par, ...)
  mp <- x$modelpar
  if(!is.null(mp)) {
    cat(paste("Derived parameters of",
              if(!is.null(mo)) mo else "model",
              "($modelpar):\n"))
    print(mp)
  }
  # Diagnostics
  printStatus(optimStatus(x$opt))
  # Starting values
  cat("Starting values of parameters:\n")
  print(x$startpar)
  # Algorithm parameters
  ct <- x$ctrl
  cat(paste("Domain of integration:",
            "[",
            signif(ct$rmin,4),
            ",",
            signif(ct$rmax,4),
            "]\n"))
  cat(paste("Exponents:",
            "p=", paste(signif(ct$p, 3), ",",  sep=""),
            "q=", signif(ct$q,3), "\n"))
  invisible(NULL)
}
              

plot.minconfit <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  do.call("plot.fv",
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

######  convergence status of 'optim' object

optimStatus <- function(x, call=NULL) {
  cgce <- x$convergence
  switch(paste(cgce),
         "0" = {
           simpleMessage(
                         paste("Converged successfully after",
                               x$counts[["function"]],
                               "iterations"),
                         call)
         },
         "1" = simpleWarning("Iteration limit maxit was reached", call),
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

############### applications (specific models) ##################


# lookup table of explicitly-known K functions and pcf
# and algorithms for computing sensible starting parameters

.Spatstat.ClusterModelInfoTable <- 
  list(
       Thomas=list(
         # Thomas process: par = (kappa, sigma2)
         modelname = "Thomas process",
         isPCP=TRUE,
         # 
         K = function(par,rvals, ...){
           if(any(par <= 0))
             return(rep(Inf, length(rvals)))
           pi*rvals^2+(1-exp(-rvals^2/(4*par[2])))/par[1]
         },
         pcf= function(par,rvals, ...){
           if(any(par <= 0))
             return(rep(Inf, length(rvals)))
           1 + exp(-rvals^2/(4 * par[2]))/(4 * pi * par[1] * par[2])
         },
         # sensible starting parameters
         selfstart = function(X) {
           kappa <- intensity(X)
           sigma2 <- 4 * mean(nndist(X))^2
           c(kappa=kappa, sigma2=sigma2)
         },
         # meaningful model parameters
         interpret = function(par, lambda) {
           kappa <- par[["kappa"]]
           sigma <- sqrt(par[["sigma2"]])
           mu <- if(is.numeric(lambda) && length(lambda) == 1)
             lambda/kappa else NA
           c(kappa=kappa, sigma=sigma, mu=mu)
         }
         ),
       # ...............................................
       MatClust=list(
         # Matern Cluster process: par = (kappa, R)
         modelname = "Matern cluster process",
         isPCP=TRUE,
         K = function(par,rvals, ..., funaux){
           if(any(par <= 0))
             return(rep(Inf, length(rvals)))
           hfun <- funaux$Hfun
           pi * rvals^2 + (1/par[1]) * hfun(rvals/(2 * par[2]))
         },
         pcf= function(par,rvals, ..., funaux){
             if(any(par <= 0))
               return(rep(Inf, length(rvals)))
             doh <- funaux$DOH
             y <- (1/(4*pi*rvals * par[1] * par[2])) * doh(rvals/(2 * par[2]))
             ifelse(rvals < .Machine$double.eps, 1, 1+y)
           },
         funaux=list(
           Hfun=function(zz) {
             ok <- (zz < 1)
             h <- numeric(length(zz))
             h[!ok] <- 1
             z <- zz[ok]
             h[ok] <- 2 + (1/pi) * (
                                    (8 * z^2 - 4) * acos(z)
                                    - 2 * asin(z)
                                    + 4 * z * sqrt((1 - z^2)^3)
                                    - 6 * z * sqrt(1 - z^2)
                                    )
             return(h)
           },
           DOH=function(zz) {
             ok <- (zz < 1)
             h <- numeric(length(zz))
             h[!ok] <- 0
             z <- zz[ok]
             h[ok] <- (16/pi) * (z * acos(z) - (z^2) * sqrt(1 - z^2))
             return(h)
           }),
         # sensible starting paramters
         selfstart = function(X) {
           kappa <- intensity(X)
           R <- 2 * mean(nndist(X)) 
           c(kappa=kappa, R=R)
         },
         # meaningful model parameters
         interpret = function(par, lambda) {
           kappa <- par[["kappa"]]
           R     <- par[["R"]]
           mu    <- if(is.numeric(lambda) && length(lambda) == 1)
             lambda/kappa else NA           
           c(kappa=kappa, R=R, mu=mu)
         }
         ),
       # ...............................................
       Cauchy=list(
         # Neyman-Scott with Cauchy clusters: par = (kappa, eta2)
         modelname = "Neyman-Scott process with Cauchy kernel",
         isPCP=TRUE,
         K = function(par,rvals, ...){
           if(any(par <= 0))
             return(rep(Inf, length(rvals)))
           pi*rvals^2 + (1 - 1/sqrt(1 + rvals^2/par[2]))/par[1]
         },
         pcf= function(par,rvals, ...){
           if(any(par <= 0))
             return(rep(Inf, length(rvals)))
           1 + ((1 + rvals^2/par[2])^(-1.5))/(2 * pi * par[2] * par[1])
         },
         selfstart = function(X) {
           kappa <- intensity(X)
           eta2 <- 4 * mean(nndist(X))^2
           c(kappa = kappa, eta2 = eta2)
         },
         # meaningful model parameters
         interpret = function(par, lambda) {
           kappa <- par[["kappa"]]
           omega <- sqrt(par[["eta2"]])/2
           mu <- if(is.numeric(lambda) && length(lambda) == 1)
             lambda/kappa else NA
           c(kappa=kappa, omega=omega, mu=mu)
         }
         ),
       # ...............................................
       VarGamma=list(
         # Neyman-Scott with VarianceGamma/Bessel clusters: par = (kappa, eta)
         modelname = "Neyman-Scott process with Variance Gamma kernel",
         isPCP=TRUE,
         K = local({
           # K function requires integration of pair correlation
           xgx <- function(x, par, nu.pcf) {
             # x * pcf(x) without check on par values
             numer <- (x/par[2])^nu.pcf * besselK(x/par[2], nu.pcf)
             denom <- 2^(nu.pcf+1) * pi * par[2]^2 * par[1] * gamma(nu.pcf + 1)
             return(x * (1 + numer/denom))
           }
           vargammaK <- function(par,rvals, ..., margs){
             # margs = list(.. nu.pcf.. ) 
             if(any(par <= 0))
               return(rep(Inf, length(rvals)))
             nu.pcf <- margs$nu.pcf
             out <- numeric(length(rvals))
             ok <- (rvals > 0)
             rvalsok <- rvals[ok]
             outok <- numeric(sum(ok))
             for (i in 1:length(rvalsok))
               outok[i] <- 2 * pi * integrate(xgx,
                                              lower=0, upper=rvalsok[i],
                                              par=par, nu.pcf=nu.pcf)$value
             out[ok] <- outok
             return(out)
           }
           vargammaK
           }), # end of 'local'
         pcf= function(par,rvals, ..., margs){
           # margs = list(..nu.pcf..)
           if(any(par <= 0))
             return(rep(Inf, length(rvals)))
           nu.pcf <- margs$nu.pcf
           sig2 <- 1 / (4 * pi * (par[2]^2) * nu.pcf * par[1])
           denom <- 2^(nu.pcf - 1) * gamma(nu.pcf)
           rr <- rvals / par[2]
           # Matern correlation function
           fr <- ifelse(rr > 0,
                        (rr^nu.pcf) * besselK(rr, nu.pcf) / denom,
                        1)
           return(1 + sig2 * fr)
         },
         parhandler = function(..., nu.ker = -1/4) {
           stopifnot(nu.ker > -1/2)
           return(list(type="Kernel",
                       model="VarGamma",
                       margs=list(nu.ker=nu.ker,
                                  nu.pcf=2*nu.ker+1)))
         },
         # sensible starting values
         selfstart = function(X) {
           kappa <- intensity(X)
           eta <- 2 * mean(nndist(X))
           c(kappa=kappa, eta=eta)
         },
         # meaningful model parameters
         interpret = function(par, lambda) {
           kappa <- par[["kappa"]]
           omega <- par[["eta"]]
           mu <- if(is.numeric(lambda) && length(lambda) == 1)
             lambda/kappa else NA
           c(kappa=kappa, omega=omega, mu=mu)
         }
         ),
       # ...............................................
       LGCP=list(
         # Log Gaussian Cox process: par = (sigma2, alpha)
         modelname = "Log-Gaussian Cox process",
         isPCP=FALSE,
         # calls Covariance() from RandomFields package
         K = function(par, rvals, ..., model, margs) {
           if(any(par <= 0))
             return(rep(Inf, length(rvals)))
           if(model == "exponential") {
             # For efficiency
             integrand <- function(r,par,...) 2*pi*r*exp(par[1]*exp(-r/par[2]))
           } else {
             # RandomFields must be loaded (this is checked by parhandler)
             integrand <- function(r,par,model,margs)
               2*pi *r *exp(Covariance(r,model=model,
                                       param=c(0.0,par[1],0.0,par[2],margs)))
           }
           nr <- length(rvals)
           th <- numeric(nr)
           if(spatstat.options("fastK.lgcp")) {
             # integrate using Simpson's rule
             fvals <- integrand(r=rvals, par=par, model=model, margs=margs)
             th[1] <- rvals[1] * fvals[1]/2
             if(nr > 1)
               for(i in 2:nr)
                 th[i] <- th[i-1] +
                   (rvals[i] - rvals[i-1]) * (fvals[i] + fvals[i-1])/2
           } else {
             # integrate using 'integrate'
             th[1] <- if(rvals[1] == 0) 0 else 
             integrate(integrand,lower=0,upper=rvals[1],
                       par=par,model=model,margs=margs)$value
             for (i in 2:length(rvals)) {
               delta <- integrate(integrand,
                                  lower=rvals[i-1],upper=rvals[i],
                                  par=par,model=model,margs=margs)
               th[i]=th[i-1]+delta$value
             }
           }
           return(th)
         },
         pcf= function(par, rvals, ..., model, margs) {
           if(any(par <= 0))
             return(rep(Inf, length(rvals)))
           if(model == "exponential") {
             # For efficiency and to avoid need for RandomFields package
             gtheo <- exp(par[1]*exp(-rvals/par[2]))
           } else {
             gtheo <- exp(Covariance(rvals,model=model,
                                     param=c(0.0,par[1],0.0,par[2],margs)))
           }
           return(gtheo)
         },
         parhandler=function(model = "exponential", ...) {
           if(!is.character(model))
             stop("Covariance function model should be specified by name")
           margs <- c(...)
           if(model != "exponential") {
             if(!(require(RandomFields) && RandomFieldsSafe()))
               stop("The package RandomFields is required")
             # check validity
             ok <- try(Covariance(0, model=model,param=c(0,1,0,1,margs)))
             if(inherits(ok, "try-error"))
               stop("Error in evaluating Covariance")
           }
           return(list(type="Covariance",model=model, margs=margs))
         },
         # sensible starting values
         selfstart = function(X) {
           alpha <- 2 * mean(nndist(X))
           c(sigma2=1, alpha=alpha)
         },
         # meaningful model parameters
         interpret = function(par, lambda) {
           sigma2 <- par[["sigma2"]]
           alpha  <- par[["alpha"]]
           mu <- if(is.numeric(lambda) && length(lambda) == 1 && lambda > 0)
             log(lambda) - sigma2/2 else NA
           c(sigma2=sigma2, alpha=alpha, mu=mu)
         }
         )
  )

spatstatClusterModelInfo <- function(name) {
  if(!is.character(name) || length(name) != 1)
    stop("Argument must be a single character string", call.=FALSE)
  nama2 <- names(.Spatstat.ClusterModelInfoTable)
  if(!(name %in% nama2))
    stop(paste(sQuote(name), "is not recognised;",
               "valid names are", commasep(sQuote(nama2))),
         call.=FALSE)
  out <- .Spatstat.ClusterModelInfoTable[[name]]
  return(out)
}

getdataname <- function(defaultvalue, ..., dataname=NULL) {
  if(!is.null(dataname)) dataname else defaultvalue
}
  
thomas.estK <- function(X, startpar=c(kappa=1,sigma2=1),
                        lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...) {

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)

  if(inherits(X, "fv")) {
    K <- X
    if(!(attr(K, "fname") %in% c("K", "K[inhom]")))
      warning("Argument X does not appear to be a K-function")
  } else if(inherits(X, "ppp")) {
    K <- Kest(X)
    dataname <- paste("Kest(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("kappa","sigma2"))

  info <- spatstatClusterModelInfo("Thomas")
  theoret <- info$K
  
  result <- mincontrast(K, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)",
                          desc="minimum contrast fit of Thomas process"),
                        explain=list(dataname=dataname,
                          fname=attr(K, "fname"),
                          modelname="Thomas process"), ...)
  # imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "sigma2")
  result$par <- par
  # infer meaningful model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="Thomas")
  return(result)
}

lgcp.estK <- function(X, startpar=c(sigma2=1,alpha=1),
                      covmodel=list(model="exponential"), 
                      lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...) {

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)
  
  if(inherits(X, "fv")) {
    K <- X
    if(!(attr(K, "fname") %in% c("K", "K[inhom]")))
      warning("Argument X does not appear to be a K-function")
  } else if(inherits(X, "ppp")) {
    K <- Kest(X)
    dataname <- paste("Kest(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("sigma2","alpha"))

  info <- spatstatClusterModelInfo("LGCP")
  
  # digest parameters of Covariance model and test validity
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
  # imbue with meaning
  par <- result$par
  names(par) <- c("sigma2", "alpha")
  result$par <- par
  result$covmodel <- cmodel
  # infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="lgcp")
  return(result)
}



matclust.estK <- function(X, startpar=c(kappa=1,R=1),
                          lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...) {

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)

  if(inherits(X, "fv")) {
    K <- X
    if(!(attr(K, "fname") %in% c("K", "K[inhom]")))
      warning("Argument X does not appear to be a K-function")
  } else if(inherits(X, "ppp")) {
    K <- Kest(X)
    dataname <- paste("Kest(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("kappa","R"))

  info <- spatstatClusterModelInfo("MatClust")
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
  # imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "R")
  result$par <- par
  # infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="MatClust")
  return(result)
}


## versions using pcf (suggested by Jan Wild)

thomas.estpcf <- function(X, startpar=c(kappa=1,sigma2=1),
                          lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...,
                          pcfargs=list()){

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)

  if(inherits(X, "fv")) {
    g <- X
    if(!(attr(g, "fname") %in% c("g", "g[inhom]")))
      warning("Argument X does not appear to be a pair correlation function")
  } else if(inherits(X, "ppp")) {
    g <- do.call("pcf.ppp", append(list(X), pcfargs))
    dataname <- paste("pcf(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("kappa","sigma2"))

  info <- spatstatClusterModelInfo("Thomas")
  theoret <- info$pcf
  
  # avoid using g(0) as it may be infinite
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
  # imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "sigma2")
  result$par <- par
  # infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="Thomas")
  return(result)
}

matclust.estpcf <- function(X, startpar=c(kappa=1,R=1),
                            lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...,
                            pcfargs=list()){

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)

  if(inherits(X, "fv")) {
    g <- X
    if(!(attr(g, "fname") %in% c("g", "g[inhom]")))
      warning("Argument X does not appear to be a pair correlation function")
  } else if(inherits(X, "ppp")) {
    g <- do.call("pcf.ppp", append(list(X), pcfargs))
    dataname <- paste("pcf(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("kappa","R"))

  info <- spatstatClusterModelInfo("MatClust")
  theoret <- info$pcf
  funaux <-  info$funaux
  
  # avoid using g(0) as it may be infinite
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
  # imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "R")
  result$par <- par
  # infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="MatClust")
  return(result)
}

lgcp.estpcf <- function(X, startpar=c(sigma2=1,alpha=1),
                      covmodel=list(model="exponential"), 
                        lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...,
                        pcfargs=list()) {
  
  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)
  
  if(inherits(X, "fv")) {
    g <- X
    if(!(attr(g, "fname") %in% c("g", "g[inhom]")))
      warning("Argument X does not appear to be a pair correlation function")
  } else if(inherits(X, "ppp")) {
    g <- do.call("pcf.ppp", append(list(X), pcfargs))
    dataname <- paste("pcf(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("sigma2","alpha"))

  info <- spatstatClusterModelInfo("LGCP")
  
  # digest parameters of Covariance model and test validity
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
  # imbue with meaning
  par <- result$par
  names(par) <- c("sigma2", "alpha")
  result$par <- par
  result$covmodel <- cmodel
  # infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="lgcp")
  return(result)
}


cauchy.estK <- function(X, startpar=c(kappa=1,eta2=1),
                        lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...) {

# omega: scale parameter of Cauchy kernel function
# eta: scale parameter of Cauchy pair correlation function
# eta = 2 * omega

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)

  if(inherits(X, "fv")) {
    K <- X
    if(!(attr(K, "fname") %in% c("K", "K[inhom]")))
      warning("Argument X does not appear to be a K-function")
  } else if(inherits(X, "ppp")) {
    K <- Kest(X)
    dataname <- paste("Kest(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("kappa","eta2"))

  info <- spatstatClusterModelInfo("Cauchy")
  theoret <- info$K

  desc <- "minimum contrast fit of Neyman-Scott process with Cauchy kernel"
  result <- mincontrast(K, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)", desc=desc),
                        explain=list(dataname=dataname,
                          fname=attr(K, "fname"),
                          modelname="Cauchy process"), ...)
  # imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "eta2")
  result$par <- par
  # infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="Cauchy")
  return(result)
}


cauchy.estpcf <- function(X, startpar=c(kappa=1,eta2=1),
                          lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...,
                          pcfargs=list()) {

# omega: scale parameter of Cauchy kernel function
# eta: scale parameter of Cauchy pair correlation function
# eta = 2 * omega

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)

  if(inherits(X, "fv")) {
    g <- X
    if(!(attr(g, "fname") %in% c("g", "g[inhom]")))
      warning("Argument X does not appear to be a pair correlation function")
  } else if(inherits(X, "ppp")) {
    g <- do.call("pcf.ppp", append(list(X), pcfargs))
    dataname <- paste("pcf(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("kappa","eta2"))

  info <- spatstatClusterModelInfo("Cauchy")
  theoret <- info$pcf

  # avoid using g(0) as it may be infinite
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
  # imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "eta2")
  result$par <- par
  # infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="Cauchy")
  return(result)
}


vargamma.estK <- function(X, startpar=c(kappa=1,eta=1), nu.ker = -1/4,
                        lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...) {

# nu.ker: smoothness parameter of Variance Gamma kernel function
# omega: scale parameter of kernel function
# nu.pcf: smoothness parameter of Variance Gamma pair correlation function
# eta: scale parameter of Variance Gamma pair correlation function
# nu.pcf = 2 * nu.ker + 1    and    eta = omega

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)

  if(inherits(X, "fv")) {
    K <- X
    if(!(attr(K, "fname") %in% c("K", "K[inhom]")))
      warning("Argument X does not appear to be a K-function")
  } else if(inherits(X, "ppp")) {
    K <- Kest(X)
    dataname <- paste("Kest(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("kappa","eta"))

  info <- spatstatClusterModelInfo("VarGamma")
  theoret <- info$K
  
  # test validity of parameter nu and digest
  ph <- info$parhandler
  cmodel <- ph(nu.ker=nu.ker)
  margs <- cmodel$margs

  desc <- "minimum contrast fit of Neyman-Scott process with Variance Gamma kernel"
  result <- mincontrast(K, theoret, startpar,
                        ctrl=list(q=q, p=p,rmin=rmin, rmax=rmax),
                        fvlab=list(label="%s[fit](r)", desc=desc),
                        explain=list(dataname=dataname,
                          fname=attr(K, "fname"),
                          modelname="Variance Gamma process"),
                        margs=margs, ...)
  # imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "eta")
  result$par <- par
  result$covmodel <- cmodel
  # infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="VarGamma")
  return(result)
}


vargamma.estpcf <- function(X, startpar=c(kappa=1,eta=1), nu.ker=-1/4, 
                          lambda=NULL, q=1/4, p=2, rmin=NULL, rmax=NULL, ...,
                          pcfargs=list()) {

# nu.ker: smoothness parameter of Variance Gamma kernel function
# omega: scale parameter of kernel function
# nu.pcf: smoothness parameter of Variance Gamma pair correlation function
# eta: scale parameter of Variance Gamma pair correlation function
# nu.pcf = 2 * nu.ker + 1    and    eta = omega

  dataname <-
    getdataname(short.deparse(substitute(X), 20), ...)

  if(inherits(X, "fv")) {
    g <- X
    if(!(attr(g, "fname") %in% c("g", "g[inhom]")))
      warning("Argument X does not appear to be a pair correlation function")
  } else if(inherits(X, "ppp")) {
    g <- do.call("pcf.ppp", append(list(X), pcfargs))
    dataname <- paste("pcf(", dataname, ")", sep="")
    if(is.null(lambda))
      lambda <- summary(X)$intensity
  } else 
    stop("Unrecognised format for argument X")

  startpar <- check.named.vector(startpar, c("kappa","eta"))

  info <- spatstatClusterModelInfo("VarGamma")
  theoret <- info$pcf

  # test validity of parameter nu and digest 
  ph <- info$parhandler
  cmodel <- ph(nu.ker=nu.ker)
  margs <- cmodel$margs
  
  # avoid using g(0) as it may be infinite
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
  # imbue with meaning
  par <- result$par
  names(par) <- c("kappa", "eta")
  result$par <- par
  result$covmodel <- cmodel
  # infer model parameters
  result$modelpar <- info$interpret(par, lambda)
  result$internal <- list(model="VarGamma")
  return(result)
}

