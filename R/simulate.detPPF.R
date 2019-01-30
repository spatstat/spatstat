##   simulate.detPPF.R
##            $Revision: 1.7 $  $Date: 2019/01/29 05:21:22 $
##
## This file contains functions to simulate DPP models.
## Two simulation functions are visible:
## - simulate.detpointprocfamily (most useful)
## - rdpp (more generic workhorse function -- actually the real workhorse is the locally defined rdppp)
##
## Furthermore the auxilliary function dppeigen is defined here.

rdpp <- local({

## Generates an empty point pattern
emptyppx <- function(W, simplify = TRUE){
  W <- as.boxx(W)
  r <- W$ranges
  d <- ncol(r)
  if(simplify){
      if(d==2)
          return(ppp(numeric(0), numeric(0), window=as.owin(W)))
      if(d==3)
          return(pp3(numeric(0), numeric(0), numeric(0), W))
  }
  rslt <- replicate(d, numeric(0), simplify=FALSE)
  names(rslt) <- paste("x",1:d,sep="")
  rslt <- as.data.frame(rslt)
  return(ppx(rslt, domain = W, coord.type= rep("spatial", d)))
}

rdppp <- function(index, basis = "fourierbasis", window = boxx(rep(list(0:1), ncol(index))),
                  reject_max = 1e4, progress = 0, debug = FALSE, given = NULL, given_max_volume = 0.5, ...){
  ## Check arguments:
  if (!(is.logical(debug)))
    stop(paste(sQuote("debug"), "must be TRUE or FALSE"))
  if (!is.numeric(reject_max)||reject_max<=1)
    stop(paste(sQuote("reject_max"), "must be a numeric greater than 1"))
  if (!is.numeric(progress)||reject_max<1)
    stop(paste(sQuote("progress"), "must be a numeric greater than or equal to 1"))
  index <- as.matrix(index)
  d <- ncol(index)
  window <- as.boxx(window)
  ranges <- window$ranges
  boxlengths <- as.numeric(ranges[2L, ] - ranges[1L, ])
  if(ncol(ranges)!=d)
    stop("The dimension differs from the number of columns in index")
  if(basis != "fourierbasis"){
    warning("Non Fourier basis probably doesn't work correctly! Fourier is
            assumed for bounds in rejection sampling.")
    userbasis <- get(basis)
    if (!(is.function(userbasis)))
      stop(paste(sQuote("basis"), "must be a function"))
    tmp <- userbasis(ranges[1,,drop=FALSE], index, window)
    if (!(is.numeric(tmp) || is.complex(tmp)))
      stop(paste("Output of", sQuote("basis"), "must be numeric or complex"))
    basis <- function(x, k, boxlengths){
      userbasis(x, k, boxx(lapply(boxlengths, function(x) list(c(0,x)))))
    }
  } else{
    basis <- fourierbasisraw
  }

  ## Number of points to simulate:
  n <- nrow(index)
  
  ## Resolve `given` for pseudo conditional simulation
  if(!is.null(given)){
    # Make sure `given` is a list of point patterns
    if(is.ppp(given) || is.pp3(given) || is.ppx(given)){
      given <- list(given)
    }
    stopifnot(all(sapply(given, function(x){ is.ppp(x) || is.pp3(x) || is.ppx(x)})))
    # Check that the window (or its boundingbox) is inside the simulation window
    Wgiven <- lapply(given, function(x) as.boxx(domain(x)))
    stopifnot(all(sapply(Wgiven,
      function(w){
        all(w$ranges[1,] >= ranges[1,]) && all(w$ranges[2,] <= ranges[2,])
        })))
    stopifnot(sum(sapply(Wgiven, volume))<given_max_volume)
    # Resolve number of given points and extract coordinates
    ngiven <- sum(sapply(given, npoints))
    stopifnot(ngiven <= n)
    if(ngiven == n) return(given)
    coordsgiven <- lapply(given, function(x) as.matrix(coords(x)))
    coordsgiven <- Reduce(rbind, coordsgiven)
  }
  
  ## Return empty point pattern if n=0:
  empty <- emptyppx(window)
  if (n==0)
    return(empty)

  ## Initialize debug info:
  if(debug){
    debugList = replicate(n, list(old=empty, accepted=empty, rejected=empty, index=index), simplify=FALSE)
  }
  
  # Matrix of coordinates:
  x <- matrix(0,n,d)
  colnames(x) <- paste("x",1:d,sep="")
  x[1,] <- runif(d,as.numeric(ranges[1,]),as.numeric(ranges[2,]))
  if(!is.null(given)){
    x[1,] <- coordsgiven[1,,drop=FALSE]
  }
  
  # Debug info:
  if(debug){
    debugList[[1]]=list(old=empty, accepted=ppx(x[1,,drop=FALSE],window,simplify=TRUE), rejected=empty, index=index, estar=rep(1/n,n))
  }
  
  if (n==1)
    return(ppx(x, window, simplify = TRUE))
  
  # First vector of basis-functions evaluated at first point:
  v <- basis(x[1,,drop=FALSE],index,boxlengths)
  ## Record normalized version in the Gram-Schmidt matrices:
  e <- v/sqrt(sum(abs(v)^2))
  estar <- Conj(e)
  if(progress>0)
    cat(paste("Simulating", n, "points:\n"))

  ## Main for loop over number of points:
  for(i in (n-1):1){
    ## Print progress:
    if(progress>0)
      progressreport(n-i, n, every=progress)
    ## Aux. variable to count number of rejection steps:
    tries <- 0
    # Debug info:
    if(debug){
      rejected <- matrix(NA,reject_max,d)
    }
    repeat{
      ## Proposed point:
      newx <- matrix(runif(d,as.numeric(ranges[1,]),as.numeric(ranges[2,])),ncol=d)
      if(!is.null(given)){
        if(i>(n-ngiven)){
          newx <- coordsgiven[n-i+1,,drop=FALSE]
        } else{
          while(any(sapply(Wgiven, function(w) inside.boxx(as.hyperframe(newx), w = w))))
            newx <- matrix(runif(d,as.numeric(ranges[1,]),as.numeric(ranges[2,])),ncol=d)
        }
      }
      ## Basis functions eval. at proposed point:
      v <- as.vector(basis(newx, index, boxlengths))
      ## Vector of projection weights (has length n-i)
      wei <- t(v)%*%estar
      if(!is.null(given) && i>(n-ngiven)){
        break
      }
      ## Accept probability:
      # tmp <- prod(ranges[2,]-ranges[1,])/n*(sum(abs(v)^2)-sum(abs(wei)^2))
      tmp <- 1-prod(ranges[2,]-ranges[1,])/n*(sum(abs(wei)^2))
      ## If proposal is accepted the loop is broken:
      if(runif(1)<as.numeric(abs(tmp))){
        break
      }
      ## If rejected, check that we have not tried too many times:
      if(tries>reject_max){
        stop(paste("Rejection sampling failed reject_max =",reject_max,"times in a row"))
      }
      ## Increase the count of rejection steps:
      tries <- tries+1
      # Debug info:
      if(debug){
        rejected[tries,] <- newx
      }
    } ## END OF REJECTION LOOP

    # Record the accepted point:
    x[n-i+1,] <- newx

    # Debug info:
    if(debug){
      if(tries==0){
        rej <- empty
      } else{
        rej <- ppx(rejected[1:tries,,drop=FALSE],window, simplify=TRUE)
      }
      debugList[[n-i+1]] = list(
                old=ppx(x[1:(n-i),,drop=FALSE],window, simplify=TRUE),
                accepted=ppx(newx,window,simplify=TRUE),
                rejected=rej, index=index, estar = estar)
    }

    ## If it is the last point exit the main loop:
    if(i==1){break}

    ## Calculate orthogonal vector for Gram-Schmidt procedure:
    # w <- v - rowSums(matrix(wei,n,n-i,byrow=TRUE)*e[,1:(n-i)])
    w <- v - colSums(t(e)*as.vector(wei))
    ## Record normalized version in the Gram-Schmidt matrices:
    enew <- w/sqrt(sum(abs(w)^2))
    e <- cbind(e, enew)
    estar <- cbind(estar,Conj(enew))
  } ## END OF MAIN FOR LOOP
  # Save points as point pattern:
  X <- ppx(x, window, simplify = TRUE)
  # Debug info:
  if(debug){
    attr(X, "dpp") <- list(debug=debugList)
  }
  if(progress>0)
    cat(" Done!\n")
  return(X)
}

rdpp <- function(eig, index, basis = "fourierbasis",
                 window = boxx(rep(list(0:1), ncol(index))), reject_max = 1e4,
                 progress = 0, debug = FALSE, ...){
  window2d <- NULL
  if (is.owin(window)) 
    window2d <- window
  sampleindex <- as.matrix(index[rbinom(nrow(index), 1, eig)==1, ])
  X <- rdppp(sampleindex, basis=basis, window=window, reject_max=reject_max, progress=progress, debug=debug, ...)
  if(!is.null(window2d))
    X <- X[window2d]
  return(X)
}

rdpp
}
)

simulate.dppm <-
simulate.detpointprocfamily <- function(object, nsim = 1, seed = NULL, ..., W = NULL,
                              trunc = .99, correction = "periodic", rbord = reach(object)
                              #  parallel = FALSE
                              ){
  
  # .... copied from simulate.lm ....
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  # ..................................

  if(inherits(object, "dppm")){
      if(is.null(W))
          W <- Window(object$X)
      object <- object$fitted
  }
  if(!inherits(object, "detpointprocfamily"))
    stop("The model to simulate must be of class detpointprocfamily")
  if(length(tmp <- object$freepar)>0)
    stop(paste("The model to simulate must be completely specified. The following parameters are unspecified:", tmp))
  if(!valid(object))
    stop("The model is invalid. Please change parameter values to get a valid model")
  if(!is.numeric(nsim)||nsim<1)
    stop(paste(sQuote("nsim"), "must be a numeric greater than or equal to 1"))
  nsim <- floor(nsim)
  dim <- dim(object)
  basis <- object$basis
  ####### BACKDOOR TO SPHERICAL CASE ########
  if(!is.null(spherefun <- object$sim_engine)){
      sphereSimEngine <- get(spherefun)
      rslt <- sphereSimEngine(object, trunc, nsim, ...)
      attr(rslt, "seed") <- RNGstate
      return(rslt)
  }
  ###########################################

  # Check stationarity and window compatibility (if 'W' and 'thin' both are present)
  statmodel <- is.null(thin <- object$thin)
  if(is.null(W)){
    if(!statmodel) W <- domain(thin)
  }
  Wowin <- if(is.owin(W)) W else NULL
  if(is.null(W)){
    W <- boxx(rep(list(0:1), dim))
  } else{
    W <- as.boxx(W, warn.owin = FALSE)
  }
  if(!statmodel){
    if(!is.subset.owin(Wowin,thin))
      stop("The window of simulation is not contained in the window of the inhomogeneous intensity.")
  }
  r <- W$ranges
  if(dim!=ncol(r))
    stop(paste("The dimension of the window:", ncol(r), "is inconsistent with the dimension of the model:", dim))
  Wscale <- as.numeric(r[2,]-r[1,])
  Wcenter <- as.numeric(colMeans(r))
  if(correction=="border"){
    if(!is.numeric(rbord)||any(rbord<0))
      stop(paste(sQuote("rbord"), "must be a non-negative numeric"))
    borderscale <- pmin((Wscale+2*rbord)/Wscale, 2)
    Wscale <- borderscale*Wscale
  }
  
  ##  lambda <- intensity(object)
  tmp <- dppeigen(object, trunc, Wscale)
  trunc <- tmp$trunc
  prec <- tmp$prec
  n <- length(tmp$eig)
  indexlist <- replicate(nsim, {x <- as.matrix(tmp$index[rbinom(n, 1, tmp$eig)==1, ]); gc(); x}, simplify = FALSE)
  rm(tmp)
  gc()
  onesim <- function(i, win=NULL){
    X <- rdpp(1, indexlist[[i]], basis = basis, window = boxx(rep(list(c(-.5,.5)), dim)), ...)
    a <- attr(X, "dpp")
    a <- c(a, list(prec = prec, trunc = trunc))
    if(correction=="border"){
      if(dim!=2)
        stop("Border correction only implemented for dimension 2 at the moment.")
      X <- X[affine.owin(as.owin(X), mat = diag(1/borderscale))]
    }
    if(is.ppp(X)){
      X <- affine(X, matrix(c(Wscale[1],0,0,Wscale[2]), 2, 2), Wcenter)
      if(!is.null(win))
        X <- X[win]
    } else{
      X <- ppx(X$data, domain = as.boxx(X$domain), coord.type = rep("spatial", dim))
      X$data <- as.hyperframe(as.data.frame(X$data)*matrix(Wscale, nrow(X$data), ncol(X$data), byrow = TRUE))
      X$domain$ranges <- X$domain$ranges*matrix(Wscale, 2, dim, byrow = TRUE) + matrix(Wcenter, 2, dim, byrow = TRUE)
      X <- ppx(X$data, X$domain, simplify = TRUE)
    }
    attr(X, "dpp") <- a
    return(X)
  }
  if(nsim==1){
    rslt <- onesim(1,win=Wowin)
    if(!statmodel)
        rslt <- rthin(rslt, P=thin)
  } else{
######## Old code for parallel simulation #########
#     if(is.logical(parallel)){
#       cl.cores <- if(parallel) NULL else 1
#     } else{
#       cl.cores <- parallel
#     }
#     rslt <- detlapply(1:nsim, onesim, cl.cores=cl.cores, win=Wowin)
###################################################
    rslt <- lapply(1:nsim, onesim, win=Wowin)
    if(!statmodel)
        rslt <- lapply(rslt, rthin, P=thin)
    names(rslt) <- paste("Simulation", 1:nsim)
    rslt <- if(dim == 2) as.solist(rslt) else as.anylist(rslt)
  }
  attr(rslt, "seed") <- RNGstate
  return(rslt)

}

dppeigen <- function(model, trunc, Wscale, stationary = FALSE){
    dim <- dim(model)
    if(stationary && dim!=2)
        stop("Stationarity can only be exploited in dimension 2 at the moment.")
    Wscale <- as.numeric(Wscale)
    check.nvector(Wscale, dim, things="dimensions")
    ## Calculate expected number of points if the intensity is a parameter
    expnum <- NULL
    lambdaname <- model$intensity
    if(!is.null(lambdaname))
        expnum <- getElement(model$fixedpar, lambdaname)*prod(Wscale)
    ## Get the maximal truncation in each dimension
    maxtrunc <- spatstat.options("dpp.maxmatrix")^(1/dim)
    ## Extract spectral density
    specden <- dppspecden(model)
    truncrange <- dppspecdenrange(model)*max(Wscale)
    
    if(trunc>=1){ ## Integer truncation fixed by user.
        if(stationary){
             ## Coordinates on axes:
            index1a <- c(rep(0,trunc),1:trunc)
            index2a <- c(1:trunc,rep(0,trunc))
            ## Coordinates of ordinary points:
            index1 <- rep(1:trunc,trunc)
            index2 <- rep(1:trunc,each=trunc)
            ## Spectral densities:
            eigo <- specden(0)
            eiga <- specden(sqrt((index1a/Wscale[1])^2+(index2a/Wscale[2])^2))
            eig <- specden(sqrt((index1/Wscale[1])^2+(index2/Wscale[2])^2))
            prec <- (eigo+2*sum(eiga)+4*sum(eig))/expnum   
        } else{
            trunc <- floor(trunc)
            index <- do.call(expand.grid, replicate(dim, seq(-trunc,trunc), simplify=FALSE))
            indexscaled <- index*matrix(1/Wscale, nrow(index), ncol(index), byrow = TRUE)
            if(model$isotropic){
                eig <- specden(sqrt(rowSums(indexscaled^2)))
            } else{
                eig <- specden(indexscaled)
            }
            prec <- sum(eig)/expnum
        }
    } else{ ## Integer truncation calculated from user-specified precision.
        if(is.null(expnum))
            stop("Cannot calculate truncation adaptively in a model without intensity parameter. Please specify trunc directly as a positive integer.")
        prec0 <- trunc
        trunc <- 1
        prec <- 0
        ## cat("truncation is being calculated adaptively. Current truncation:\n")
        while(prec<=prec0 && (2*trunc)<=maxtrunc && trunc<=truncrange){
            trunc <- 2*trunc
            if(stationary){
                ## Coordinates on axes:
                index1a <- c(rep(0,trunc),1:trunc)
                index2a <- c(1:trunc,rep(0,trunc))
                ## Coordinates of ordinary points:
                index1 <- rep(1:trunc,trunc)
                index2 <- rep(1:trunc,each=trunc)
                ## Spectral densities:
                eigo <- specden(0)
                eiga <- specden(sqrt((index1a/Wscale[1])^2+(index2a/Wscale[2])^2))
                eig <- specden(sqrt((index1/Wscale[1])^2+(index2/Wscale[2])^2))
                prec <- (eigo+2*sum(eiga)+4*sum(eig))/expnum
            } else{
                index <- do.call(expand.grid, replicate(dim, seq(-trunc,trunc), simplify=FALSE))
                indexscaled <- index*matrix(1/Wscale, nrow(index), ncol(index), byrow = TRUE)
                if(model$isotropic){
                    eig <- specden(sqrt(rowSums(indexscaled^2)))
                } else{
                    eig <- specden(indexscaled)
                }
                prec <- sum(eig)/expnum
            }
        }
        ## cat("\n")
        if(prec<prec0){
            warning(paste0("Adaptive truncation stopped at ", trunc, ". The precision is only ", prec))
        }
    }
    if(stationary){
        rslt <- list(eigo=eigo, eiga=eiga, eig=eig, index1a=index1a, index2a=index2a)
    } else{
        rslt <- list(eig=eig, index=index)
    }
    return(c(rslt, list(prec=prec, trunc=trunc)))
}
