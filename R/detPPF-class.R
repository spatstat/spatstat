## support for class 'detpointprocfamily'

print.detpointprocfamily <- function(x, ...){
  splat(x$name, "determinantal point process model",
        ifelse(is.numeric(x$dim), paste("in dimension", x$dim), ""))
  #' Not used:
  #'  parnames <- names(x$par)
  anyfixed <- length(x$fixedpar)>0
  if(anyfixed){
      fixedlambda <- NULL
      if(!is.null(x$intensity) && is.element(x$intensity, names(x$fixedpar))){
          lambda <- signif(x$fixedpar[[x$intensity]], 4)
          x$fixedpar <- x$fixedpar[names(x$fixedpar)!=x$intensity]
          fixedlambda <- paste(x$intensity, ifelse(is.null(x$thin), paste("=", lambda), "= an image"))
      }
      if(length(x$fixedpar)>0){
          fixedparstring <- paste(names(x$fixedpar), signif(unlist(x$fixed),4), sep = " = ", collapse = ", ")
          fixedparstring <- paste(fixedlambda, fixedparstring, sep=", ")
      } else{
          fixedparstring <- fixedlambda
      }
  }
  ## Partially specified model:
  if(length(x$freepar)>0){
    splat("The model is only partially specified.")
    splat("The following parameters are free (e.g. to be estimated by dppm):")
    cat(x$freepar, sep = ", ")
    cat("\n")
    if(anyfixed){
        cat("The fixed parameters are: ")
        cat(fixedparstring, sep = ", ")
    } else{
        splat("There are no fixed parameters.")
    }        
  } else{
    cat("The parameters are: ")
    cat(fixedparstring, sep = ", ")
  }
  cat("\n")
  if(!is.null(x$intensity)){
    splat("The parameter", x$intensity,
          "specifies the intensity of the process.")
  }
  if(is.character(x$dim)){
    splat("The parameter", x$dim,
          "specifies the dimension of the state space.")
  }
  invisible(NULL)
}

reach.detpointprocfamily <- function(x, ...){
    model <- x
    fun <- model$range
    nam <- names(formals(fun))
    do.call(model$range, c(model$fixedpar[is.element(names(model$fixedpar),nam)], list(...)))
}

dppparbounds <- function(model, name, ...){
    if(inherits(model, "dppm"))
        model <- model$fitted
    if(!inherits(model, "detpointprocfamily"))
        stop("input model must be of class detpointprocfamily or dppm")
    fun <- model$parbounds
    nam <- names(formals(fun))
    if(missing(name))
        name <- nam[!is.element(nam, c("name", model$dim))]
    rslt <- matrix(0,length(name), 2, dimnames = list(name, c("lower", "upper")))
    for(nn in name){
        tmp <- try(do.call(fun, c(model$fixedpar[is.element(names(model$fixedpar),nam)], list(...), list(name=nn))), silent=TRUE)
        if(class(tmp)=="try-error"){
            rslt[nn,] <- c(NA, NA)
        }else{
            rslt[nn,] <- tmp
        }
    }
    rslt
}

valid.detpointprocfamily <- function(object, ...){
  if(length(object$freepar)>0)
      return(NA)
  ## If there is no function for checking validity we always return TRUE:
  if(is.null(object$valid))
    return(TRUE)
  do.call(object$valid, object$fixedpar)
}

dppspecdenrange <- function(model){
  ## If there is no function for checking finite range of spectral density we always return Inf:
  fun <- model$specdenrange
  if(is.null(fun))
    return(Inf)
  xx <- try(fun(model), silent = TRUE)
  ifelse(class(xx)=="try-error", Inf, xx)
}

dppspecden <- function(model){
  fun <- model$specden
  if(is.null(fun))
    stop("Spectral density unknown for this model!")
  if(length(model$freepar)>0)
    stop("Cannot extract the spectral density of a partially specified model. Please supply all parameters.")
  specden <- function(x, ...){
    allargs <- c(list(x), model$fixedpar, list(...))
    do.call(fun, allargs)
  }
  return(specden)
}

dppkernel <- function(model, ...){
  if(inherits(model, "dppm"))
    model <- model$fitted
  fun <- model$kernel
  if(is.null(fun))
    return(dppapproxkernel(model, ...))
  if(length(model$freepar)>0)
    stop("Cannot extract the kernel of a partially specified model. Please supply all parameters.")
  firstarg <- names(formals(fun))[1L]
  kernel <- function(x){
    allargs <- c(structure(list(x), .Names=firstarg), model$fixedpar)
    do.call(fun, allargs)
  }
  return(kernel)
}

dppapproxkernel <- function(model, trunc = .99, W = NULL){
    if(inherits(model, "dppm")){
        W <- model$window
        model <- model$fitted
    }
    ####### BACKDOOR TO SPHERICAL CASE ########
    if(!is.null(spherefun <- model$approxkernelfun)){
        spherefun <- get(spherefun)
        rslt <- spherefun(model, trunc)
        return(rslt)
    }
    ###########################################
    d <- dim(model)
    if(is.null(W))
      W <- boxx(replicate(d, c(-.5,.5), simplify=FALSE))
    W <- as.boxx(W)
    if(d!=ncol(W$ranges))
        stop(paste("The dimension of the window:", ncol(W$ranges), "is inconsistent with the dimension of the model:", d))
    Wscale <- as.numeric(W$ranges[2L,]-W$ranges[1L,])
    tmp <- dppeigen(model, trunc, Wscale, stationary=FALSE)
    index <- tmp$index
    eig <- tmp$eig
    prec <- tmp$prec
    trunc <- tmp$trunc
    rm(tmp)
    f <- function(r){
        x <- matrix(0, nrow=length(r), ncol=d)
        x[,1L] <- r
        basis <- fourierbasis(x, index, win = W)
        approx <- matrix(eig, nrow=length(eig), ncol=length(r)) * basis
        return(Re(colSums(approx)))
    }
    attr(f, "dpp") <- list(prec = prec, trunc = trunc)
    return(f)
}

pcfmodel.detpointprocfamily <- function(model, ...){
  kernel <- dppkernel(model, ...)
  f <- function(x){
    1 - (kernel(x)/kernel(0))^2
  }
  return(f)
}

dppapproxpcf <- function(model, trunc = .99, W = NULL){
  kernel <- dppapproxkernel(model, trunc = trunc, W = W)
  f <- function(x){
    1 - (kernel(x)/kernel(0))^2
  }
  attr(f, "dpp") <- attr(kernel, "dpp")
  return(f)
}

Kmodel.detpointprocfamily <- function(model, ...){
  if(length(model$freepar)>0)
    stop("Cannot extract the K function of a partially specified model. Please supply all parameters.")
  fun <- model$Kfun
  if(!is.null(fun)){
      firstarg <- names(formals(fun))[1L]
      Kfun <- function(r){
          allargs <- c(structure(list(r), .Names=firstarg), model$fixedpar)
          do.call(fun, allargs)
      }
  } else{
      pcf <- pcfmodel(model, ...)
      intfun <- function(xx){
          2*pi*xx*pcf(xx)
      }
      Kfun <- function(r){
          r <- sort(r)
          if(r[1L]<0)
              stop("Negative values not allowed in K function!")
          r <- c(0,r)
          int <- unlist(lapply(2:length(r), function(i) integrate(intfun, r[i-1L], r[i], subdivisions=10)$value))
          return(cumsum(int))
      }
  }
  return(Kfun)
}

update.detpointprocfamily <- function(object, ...){
     newpar <- list(...)
     if(length(newpar)==1L && is.list(newpar[[1L]]) && !is.im(newpar[[1L]]))
         newpar <- newpar[[1L]]
     nam <- names(newpar)
     if(length(newpar)>0&&is.null(nam))
         stop(paste("Named arguments are required. Please supply parameter values in a", sQuote("tag=value"), "form"))
     oldpar <- object$fixedpar[!is.element(names(object$fixedpar), nam)]
     thin <- object$thin
     object <- do.call(object$caller, c(newpar,oldpar))
     if(is.null(object$thin))
        object$thin <- thin
     return(object)
}

is.stationary.detpointprocfamily <- function(x){
    if(is.null(x$intensity))
        return(FALSE)
    lambda <- getElement(x$fixedpar, x$intensity)
    if(!is.null(lambda)&&is.numeric(lambda)&&is.null(x$thin))
        return(TRUE)
    return(FALSE)
    
}

intensity.detpointprocfamily <- function(X, ...){
    lambda <- NULL
    if(!is.null(X$intensity))
        lambda <- getElement(X$fixedpar, X$intensity)
    if(!is.null(lambda)){
      if(!is.null(X$thin))
        lambda <- lambda*X$thin
      return(lambda)
    }
    return(NA)
}

parameters.dppm <- parameters.detpointprocfamily <- function(model, ...){
    if(inherits(model, "dppm"))
        model <- model$fitted
    c(model$fixed, structure(rep(NA,length(model$freepar)), .Names = model$freepar))
}

dim.detpointprocfamily <- function(x){
    if(is.numeric(d <- x$dim)){
        return(d)
    } else{
        return(getElement(x$fixedpar, d))
    }
}
