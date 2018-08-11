##    detpointprocfamilyfun.R
##
##    $Revision: 1.5 $   $Date: 2015/10/19 02:27:17 $
##
## This file contains the function `detpointprocfamilyfun'
## to define new DPP model family functions
## and a print method for class `detpointprocfamilyfun'
## as well as the currently defined 
## - dppBessel
## - dppCauchy
## - dppGauss
## - dppMatern
## - dppPowerExp

detpointprocfamilyfun <- local({

names_formals <- function(f, dots = FALSE){
    nam <- names(formals(f))
    if(!dots) nam <- nam[nam!="..."]
    return(nam)
}

detpointprocfamilyfun <-
  function(kernel=NULL, specden=NULL, basis="fourierbasis",
           convkernel=NULL, Kfun=NULL, valid=NULL,
           intensity=NULL, dim=2, name="User-defined",
           isotropic=TRUE, range=NULL, parbounds=NULL,
           specdenrange=NULL, startpar=NULL, ...)
{
  ## Check which functions are given, check them for sanity and
  ## extract argument names and other stuff
  given <- NULL
  if(!is.null(kernel)){
    if(!is.function(kernel))
      stop("If kernel is given it must be a function.")
    given <- "kernel"
    kernelnames <- names_formals(kernel)
    if(length(kernelnames)<1L)
      stop("kernel function must have at least one argument")
    kernelnames <- kernelnames[-1L]
  }
  if(!is.null(specden)){
    if(!is.function(specden))
      stop("If specden is given it must be a function.")
    given <- c(given, "specden")
    specdennames <- names_formals(specden)
    if(length(specdennames)<1L)
      stop("specden function must have at least one argument")
    specdennames <- specdennames[-1L]
  }
  if(is.null(given))
    stop("At least one of kernel or specden must be provided.")
  if(length(given)==2){
    if(!setequal(kernelnames,specdennames))
      stop("argument names of kernel and specden must match.")
  }
  if(is.element("kernel",given)){
    parnames <- kernelnames
  } else{
    parnames <- specdennames
  }
  if(!is.null(convkernel)){
    given <- c(given,"convkernel")
    if(!is.function(convkernel)||length(formals(convkernel))<2)
      stop("If convkernel is given it must be a function with at least two arguments.")
    if(!setequal(parnames,names_formals(convkernel)[-(1:2)]))
      stop("argument names of convkernel must match argument names of kernel and/or specden.")
  }
  if(!is.null(Kfun)){
    given <- c(given,"Kfun")
    if(!is.function(Kfun)||length(formals(Kfun))<1L)
      stop("If Kfun is given it must be a function with at least one arguments.")
    if(!setequal(parnames,names_formals(Kfun)[-1L]))
      stop("argument names of Kfun must match argument names of kernel and/or specden.")
  }
  if(!is.null(valid)){
    if(!(is.function(valid)&&setequal(parnames,names_formals(valid))))
      stop("argument names of valid must match argument names of kernel and/or specden.")
  } else{
    warning("No function for checking parameter validity provided. ANY numerical value for the parameters will be accepted.")
  }
  if(!is.null(intensity)&&!(is.character(intensity)&&length(intensity)==1L&&is.element(intensity, parnames)))
    stop("argument intensity must be NULL or have length one, be of class character and match a parameter name")
      
  if(!(is.character(dim)|is.numeric(dim))|length(dim)!=1L)
    stop("argument dim must have length one and be of class character or numeric")
  if(is.character(dim)){
    if(!is.element(dim, parnames))
      stop("When dim is a character it must agree with one of the parameter names of the model")
  } else{
    dim <- round(dim)
    if(dim<1L)
      stop("When dim is a numeric it must be a positive integer")
  }

  ## Catch extra unknown args (will be appended to output object below).
  dots <- list(...)

  ## Create output object.
  out <- function(...){
    caller <- match.call()[[1L]]
    caller <- eval(substitute(caller), parent.frame())
    fixedpar <- list(...)
    nam <- names(fixedpar)
    if(length(fixedpar)>0&&is.null(nam))
      stop(paste("Named arguments are required. Please supply parameter values in a", sQuote("tag=value"), "form"))
    match <- is.element(nam, parnames)
    if(sum(!match)>0)
      warning(paste("Not all supplied argument(s) make sense. Valid arguments are: ",
                    paste(parnames, collapse = ", "),
                    ". The following supplied argument(s) will be ignored: ",
                    paste(nam[!match], collapse = ", "),
                    sep = ""))
    fixedpar <- fixedpar[match]
    
    ## Code to always fix the dimension to a numeric when calling the function #######
    if(is.character(dim) && !is.element(dim,names(fixedpar))){
         dimpar <- structure(list(2), .Names=dim)
         fixedpar <- c(fixedpar, dimpar)
    }
    
    ## Detect inhomogeneous intensity (an image), and replace by max and an image for thinning
    thin <- NULL
    if(!is.null(intensity)){
       lambda <- getElement(fixedpar, intensity)
       if(is.im(lambda)){
         lambdamax <- max(lambda)
         thin <- lambda/lambdamax
         fixedpar[[intensity]] <- lambdamax
       }
    }
      
    obj <- list(fixedpar = fixedpar,
                freepar = parnames[!is.element(parnames,names(fixedpar))],
                kernel = kernel,
                specden = specden,
                convkernel = convkernel,
                intensity = intensity,
                thin = thin,
                dim = dim,
                name = name,
                range = range,
                valid = valid,
                parbounds = parbounds,
                specdenrange = specdenrange,
                startpar = startpar,
                isotropic = isotropic,
                caller = caller,
                basis = basis
                )
    obj <- append(obj, dots)
    class(obj) <- "detpointprocfamily"
    return(obj)
  }
  class(out) <- c("detpointprocfamilyfun",
                  "pointprocfamilyfun",
                  class(out))
  attr(out, "parnames") <- parnames
  attr(out, "name") <- name
  return(out)
}

detpointprocfamilyfun
}
)

print.detpointprocfamilyfun <- function(x, ...){
  cat(paste(attr(x, "name"), "determinantal point process model family\n"))
  cat("The parameters of the family are:\n")
  cat(attr(x, "parnames"), sep = ", ")
  cat("\n")
  invisible(NULL)
}

dppBessel <- detpointprocfamilyfun(
  name="Bessel",
  kernel=function(x, lambda, alpha, sigma, d){
    a <- 0.5*(sigma+d)
    y <- abs(x/alpha)
    # Kernel: lambda*2^a*gamma(a+1)*besselJ(2*y*sqrt(a),a) / (2*y*sqrt(a))^a
    logrslt <- log(lambda) + a*log(2) + lgamma(a+1) - a*log(2*y*sqrt(a))
    rslt <- exp(logrslt) * besselJ(2*y*sqrt(a), a)
    rslt[x==0] <- lambda
    return(rslt)
  },
  specden=function(x, lambda, alpha, sigma, d){
    a <- sigma+d
    # specden: lambda*(2*pi)^(d/2)*alpha^d*gamma(0.5*a+1)/a^(d/2)/gamma(sigma/2+1)*(1-2*pi^2*alpha^2*x^2/a)^(sigma/2)
    logrslt <- log(lambda) + (d/2)*log(2*pi) + d*log(alpha) + lgamma(0.5*a+1)
    logrslt <- logrslt - (d/2)*log(a) - lgamma(sigma/2+1)
    tmp <- 1-2*pi^2*alpha^2*x^2/a
    warnopt <- options(warn=-1)
    logrslt <- logrslt + ifelse(tmp<0, -Inf, (sigma/2)*log(tmp))
    options(warnopt)
    return(exp(logrslt))
  },
  specdenrange=function(model){
    p <- model$fixedpar
    sqrt((p$sigma+p$d)/(2*pi^2*p$alpha^2))
  },
  valid=function(lambda, alpha, sigma, d){
    a <- sigma+d
    OK <- lambda>0 && alpha>0 && d>=1 && sigma>=0
    if(!OK)
      return(FALSE)
    ## Upper bound for alpha (using log-scale)
    lognum <- log(a^(0.5*d)) + lgamma(0.5*sigma+1)
    logdenom <- log( lambda*(2*pi^(0.5*d))) + lgamma(0.5*a+1)
    logalphamax <- (1/d) * (lognum - logdenom)
    return(OK && log(alpha) <= logalphamax)
  },
  isotropic=TRUE,
  intensity="lambda",
  dim="d",
  parbounds=function(name, lambda, alpha, sigma, d){
    lognum <- log((sigma+d)^(0.5*d)) + lgamma(0.5*sigma+1)
    logdenom <- log(2*pi^(0.5*d)) + lgamma(0.5*(sigma+d)+1)
    switch(name,
           lambda = c(0, exp(lognum - log( alpha^d) - logdenom)) ,
           alpha = c(0, exp((1/d) * (lognum - log(lambda) - logdenom))),
           sigma = c(0, switch(as.character(d), "2"=Inf, NA)),
           stop("Parameter name misspecified")
    )
  },
  startpar=function(model, X){
    rslt <- NULL
    if("d" %in% model$freepar){
      model <- update(model, d=spatdim(X))
    }
    if("lambda" %in% model$freepar){
      lambda <- intensity(X)
      while(!is.na(OK <- valid(model <- update(model, lambda=lambda)))&&!OK)
        lambda <- lambda/2
      rslt <- c(rslt, "lambda" = lambda)
    }
    if("sigma" %in% model$freepar){
      sigma <- 2
      while(!is.na(OK <- valid(model <- update(model, sigma=sigma)))&&!OK)
        sigma <- sigma/2
      rslt <- c(rslt, "sigma" = sigma)
    }
    if("alpha" %in% model$freepar){
      alpha <- .8*dppparbounds(model, "alpha")[2L]
      while(!is.na(OK <- valid(model <- update(model, alpha=alpha)))&&!OK){
        alpha <- alpha/2
      }
      rslt <- c(rslt, "alpha" = alpha)
    }
    return(rslt)
  }
)

dppCauchy <- detpointprocfamilyfun(
  name="Cauchy",
  kernel=function(x, lambda, alpha, nu, d){
    rslt <- lambda * (1+(x/alpha)^2)^(-nu-d/2)
    rslt[x==0] <- lambda
    return(rslt)
  },
  specden=function(x, lambda, alpha, nu, d){
    y <- 2*x*alpha*pi
    rslt <- lambda * y^nu * besselK(x = y, nu = nu) * (sqrt(pi)*alpha)^d * exp((1-nu)*log(2) - lgamma(nu+d/2))
    rslt[x==0] <- lambda * exp(lgamma(nu) - lgamma(nu+d/2)) * (sqrt(pi)*alpha)^d
    return(rslt)
  },
  Kfun = function(x, lambda, alpha, nu, d){
    rslt <- pi*x^2 - pi*alpha^2/(2*nu+1) * (1 - (alpha^2/(alpha^2+x^2))^(2*nu+1))
    rslt[rslt<0] <- 0
    return(rslt)
  },
  valid=function(lambda, alpha, nu, d){
    ## Note the upper bound on nu for numerical stability!
    lambda>0 && alpha>0 && nu>0 && nu<=50 && d>=1 && lambda <= gamma(nu+d/2)/(gamma(nu)*(sqrt(pi)*alpha)^d)
  },
  isotropic=TRUE,
  intensity="lambda",
  dim="d",
  range=function(alpha, nu, d, bound = .99){
    if(missing(alpha))
      stop("The parameter alpha is missing.")
    if(missing(nu))
      stop("The parameter nu is missing.")
    if(missing(d))
      stop("The parameter d (giving the dimension) is missing.")
    if(!(is.numeric(bound)&&bound>0&&bound<1))
      stop("Argument bound must be a numeric between 0 and 1.")
    return(alpha * sqrt((1-bound)^(-1/(2*nu+d))-1))
  },
  parbounds=function(name, lambda, alpha, nu, d){
    switch(name,
           lambda = c(0, gamma(nu+d/2)/(gamma(nu)*(sqrt(pi)*alpha)^d)),
           alpha = c(0, (exp(lgamma(nu+d/2)-lgamma(nu))/lambda)^(1/d)/sqrt(pi)),
           ## nu bound only implemented for d = 2.
           nu = c(switch(as.character(d), "2"=pi*lambda*alpha^2, NA), Inf),
           stop("Parameter name misspecified")
    )
  },
  startpar=function(model, X){
    rslt <- NULL
    if("lambda" %in% model$freepar){
      lambda <- intensity(X)
      while(!is.na(OK <- valid(model <- update(model, lambda=lambda)))&&!OK)
        lambda <- lambda/2
      rslt <- c(rslt, "lambda" = lambda)
    }
    if("nu" %in% model$freepar){
      nu <- 2
      while(!is.na(OK <- valid(model <- update(model, nu=nu)))&&!OK)
        nu <- nu/2
      rslt <- c(rslt, "nu" = nu)
    }
    if("alpha" %in% model$freepar){
      alpha <- .8*dppparbounds(model, "alpha")[2L]
      while(!is.na(OK <- valid(model <- update(model, alpha=alpha)))&&!OK){
        alpha <- alpha/2
      }
      rslt <- c(rslt, "alpha" = alpha)
    }
    return(rslt)
  }
)

dppGauss <- detpointprocfamilyfun(
  name="Gaussian",
  kernel=function(x, lambda, alpha, d){
    rslt <- lambda*exp(-(x/alpha)^2)
    return(rslt)
  },
  specden=function(x, lambda, alpha, d){
    lambda * (sqrt(pi)*alpha)^d * exp(-(x*alpha*pi)^2)
  },
  convkernel=function(x, k, lambda, alpha, d){
    logres <- k*log(lambda*pi*alpha^2) - log(pi*k*alpha^2) - x^2/(k*alpha^2)
    return(exp(logres))
  },
  Kfun = function(x, lambda, alpha, d){
    pi*x^2 - pi*alpha^2/2*(1-exp(-2*x^2/alpha^2))
  },
  valid=function(lambda, alpha, d){
    lambda>0 && alpha>0 && d>=1 && lambda <= (sqrt(pi)*alpha)^(-d)
  },
  isotropic=TRUE,
  intensity="lambda",
  dim="d",
  range=function(alpha, bound = .99){
    if(missing(alpha))
      stop("The parameter alpha is missing.")
    if(!(is.numeric(bound)&&bound>0&&bound<1))
      stop("Argument bound must be a numeric between 0 and 1.")
    return(alpha*sqrt(-log(sqrt(1-bound))))
  },
  parbounds=function(name, lambda, alpha, d){
    switch(name,
           lambda = c(0, (sqrt(pi)*alpha)^(-d)),
           alpha = c(0, lambda^(-1/d)/sqrt(pi)),
           stop("Parameter name misspecified")
    )
  },
  startpar=function(model, X){
    rslt <- NULL
    if("lambda" %in% model$freepar){
      lambda <- intensity(X)
      rslt <- c(rslt, "lambda" = lambda)
      model <- update(model, lambda=lambda)
    }
    if("alpha" %in% model$freepar){
      alpha <- .8*dppparbounds(model, "alpha")[2L]
      rslt <- c(rslt, "alpha" = alpha)
    }
    return(rslt)
  }
)

dppMatern <- detpointprocfamilyfun(
  name="Whittle-Matern",
  kernel=function(x, lambda, alpha, nu, d){
    rslt <- lambda*2^(1-nu) / gamma(nu) * ((x/alpha)^nu) * besselK(x = x/alpha, nu = nu)
    rslt[x==0] <- lambda
    return(rslt)
  },
  specden=function(x, lambda, alpha, nu, d){
    lambda * exp(lgamma(nu+d/2) - lgamma(nu)) * (2*sqrt(pi)*alpha)^d * (1+(2*x*alpha*pi)^2)^(-nu-d/2)
  },
  convkernel=function(x, k, lambda, alpha, nu, d){
    nu2 <- k*(nu+d/2)-d/2
    logres <- (nu2)*log(x/alpha) + log(besselK(x = x/alpha, nu = nu2, expon.scaled = TRUE)) - x/alpha
    logres[x == 0] <- (nu2-1)*log(2) + lgamma(nu2)
    logres <- logres + k*log(lambda) + k*(lgamma(nu+d/2)-lgamma(nu)) + (d*k-d+1-nu2)*log(2) + d*(k-1)*log(sqrt(pi)*alpha) - lgamma(nu2+d/2)
    index <- which(logres == Inf)
    logres[index] <- -Inf
    return(exp(logres))
  },
  valid=function(lambda, alpha, nu, d){
    ## Note the upper bound on nu for numerical stability!
    lambda>0 && alpha>0 && nu>0 && nu<=50 && d>=1 && lambda <= gamma(nu)/(gamma(nu+d/2)*(2*sqrt(pi)*alpha)^d)
  },
  isotropic=TRUE,
  intensity="lambda",
  dim="d",
  range=function(alpha, nu, d, bound = .99, exact = FALSE){
    if(missing(alpha))
      stop("The parameter alpha is missing.")
    if(missing(nu))
      stop("The parameter nu is missing.")
    if(missing(d))
      stop("The parameter d (giving the dimension) is missing.")
    if(!is.logical(exact))
      stop("Argument exact must be a logical.")
    if(!exact&&d==2)
      return(alpha * sqrt(8*nu)) ## range suggested by Haavard Rue et al.
    if(!(is.numeric(bound)&&bound>0&&bound<1))
      stop("Argument bound must be a numeric between 0 and 1.")
    fun <- function(x) sqrt(1-bound)-2^(1-nu) / gamma(nu) * ((x/alpha)^nu) * besselK(x = x/alpha, nu = nu)
    return(uniroot(fun, c(sqrt(.Machine$double.eps),1e3*alpha*sqrt(nu)))$root)
  },
  parbounds=function(name, lambda, alpha, nu, d){
    switch(name,
           lambda = c(0, gamma(nu)/(gamma(nu+d/2)*(2*sqrt(pi)*alpha)^d)),
           alpha = c(0, (exp(lgamma(nu)-lgamma(nu+d/2))/lambda)^(1/d)/2/sqrt(pi)),
           ## nu bound only implemented for d = 2 and d = 4.
           nu = c(0, switch(as.character(d), "2"=1/(4*pi*lambda*alpha^2), "4"=sqrt(1/4+1/(lambda*16*pi*pi*alpha^4))-1/2, NA)),
           stop("Parameter name misspecified")
    )
  },
  startpar=function(model, X){
    rslt <- NULL
    if("lambda" %in% model$freepar){
      lambda <- intensity(X)
      while(!is.na(OK <- valid(model <- update(model, lambda=lambda)))&&!OK)
        lambda <- lambda/2
      rslt <- c(rslt, "lambda" = lambda)
    }
    if("nu" %in% model$freepar){
      nu <- 2
      while(!is.na(OK <- valid(model <- update(model, nu=nu)))&&!OK)
        nu <- nu/2
      rslt <- c(rslt, "nu" = nu)
    }
    if("alpha" %in% model$freepar){
      alpha <- .8*dppparbounds(model, "alpha")[2L]
      while(!is.na(OK <- valid(model <- update(model, alpha=alpha)))&&!OK){
        alpha <- alpha/2
      }
      rslt <- c(rslt, "alpha" = alpha)
    }
    return(rslt)
  }
)

dppPowerExp <- detpointprocfamilyfun(
  name="Power Exponential Spectral",
  specden=function(x, lambda, alpha, nu, d){
    lambda * gamma(d/2+1) * alpha^d / (pi^(d/2)*gamma(d/nu+1)) * exp(-(alpha*x)^nu)
  },
  valid=function(lambda, alpha, nu, d){
    ## Note the upper bound on nu for numerical stability!
    lambda>0 && alpha>0 && nu>0 && nu<=20 && d>=1 && lambda <= pi^(d/2)*gamma(d/nu+1) / (gamma(1+d/2)*alpha^d)
  },
  isotropic=TRUE,
  intensity="lambda",
  dim="d",
  parbounds=function(name, lambda, alpha, nu, d){
    switch(name,
           lambda = c(0, pi^(d/2)*gamma(d/nu+1) / (gamma(d/2+1)*alpha^d)),
           alpha = c(0, (pi^(d/2)*gamma(d/nu+1) / (lambda * gamma(d/2+1)))^(1/d)),
           nu = c(NA, NA),
           stop("Parameter name misspecified")
    )
  },
  startpar=function(model, X){
    rslt <- NULL
    if("lambda" %in% model$freepar){
      lambda <- intensity(X)
      while(!is.na(OK <- valid(model <- update(model, lambda=lambda)))&&!OK)
        lambda <- lambda/2
      rslt <- c(rslt, "lambda" = lambda)
    }
    if("nu" %in% model$freepar){
      nu <- 2
      while(!is.na(OK <- valid(model <- update(model, nu=nu)))&&!OK)
        nu <- nu/2
      rslt <- c(rslt, "nu" = nu)
    }
    if("alpha" %in% model$freepar){
      alpha <- .8*dppparbounds(model, "alpha")[2L]
      while(!is.na(OK <- valid(model <- update(model, alpha=alpha)))&&!OK){
        alpha <- alpha/2
      }
      rslt <- c(rslt, "alpha" = alpha)
    }
    return(rslt)
  }
)
