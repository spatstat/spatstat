##
##   randomNS.R
##
##   simulating from Neyman-Scott processes
##
##   $Revision: 1.23 $  $Date: 2015/10/21 09:06:57 $
##
##    Original code for rCauchy and rVarGamma by Abdollah Jalilian
##    Other code and modifications by Adrian Baddeley
##    Bug fixes by Abdollah, Adrian, and Rolf Turner

rNeymanScott <- 
  function(kappa, expand, rcluster, win = owin(c(0,1),c(0,1)), ...,
           lmax=NULL, nsim=1, drop=TRUE, nonempty=TRUE, saveparents=TRUE)
{
  ## Generic Neyman-Scott process
  ## Implementation for bounded cluster radius
  ##

  ## Catch old argument name rmax for expand
  if(missing(expand) && !is.null(rmax <- list(...)$rmax))
    expand <- rmax
    
  ## 'rcluster' may be
  ##
  ##     (1) a function(x,y, ...) that takes the coordinates
  ##         (x,y) of the parent point and generates a list(x,y) of offspring
  ##

  if(is.function(rcluster))
    return(rPoissonCluster(kappa, expand, rcluster, win, ...,
                           lmax=lmax, nsim=nsim, drop=drop,
                           saveparents=saveparents))

  ##     (2) a list(mu, f) where mu is a numeric value, function, or pixel image
  ##         and f is a function(n, ...) generating n i.i.d. offspring at 0,0
  
  if(!(is.list(rcluster) && length(rcluster) == 2))
    stop("rcluster should be either a function, or a list of two elements")
  win <- as.owin(win)
  mu <- rcluster[[1]]
  rdisplace <- rcluster[[2]]
  if(is.numeric(mu)) {
    ## homogeneous
    if(!(length(mu) == 1 && mu >= 0))
      stop("rcluster[[1]] should be a single nonnegative number")
    mumax <- mu
  } else if (is.im(mu) || is.function(mu)) {
      ## inhomogeneous
    if(is.function(mu)) mu <- as.im(mu, W=win, ..., strict=TRUE)
    mumax <- max(mu)
  } else stop("rcluster[[1]] should be a number, a function or a pixel image")  
  if(!is.function(rdisplace))
    stop("rcluster[[2]] should be a function")

  ## Generate parents in dilated window
  frame <- boundingbox(win)
  dilated <- grow.rectangle(frame, expand)
  if(is.im(kappa) && !is.subset.owin(dilated, as.owin(kappa)))
    stop(paste("The window in which the image",
               sQuote("kappa"),
               "is defined\n",
               "is not large enough to contain the dilation of the window",
               sQuote("win")))
  if(nonempty) {
    if(is.function(kappa)) {
      kappa <- as.im(kappa, W=dilated, ..., strict=TRUE)
      lmax <- NULL
    }
    ## intensity of parents with at least one offspring point
    kappa <- kappa * (1 - exp(-mumax))
  }
  ## generate
  parentlist <- rpoispp(kappa, lmax=lmax, win=dilated, nsim=nsim,
                        drop=FALSE, warnwin=FALSE)

  resultlist <- vector(mode="list", length=nsim)
  for(i in 1:nsim) {
    parents <- parentlist[[i]]
    
    np <- npoints(parents)
    ## generate cluster sizes
    if(np == 0) {
      ## no parents - empty pattern
      result <- ppp(numeric(0), numeric(0), window=win)
      parentid <- integer(0)
    } else {
      if(!nonempty) {
        ## cluster sizes are Poisson
        csize <- rpois(np, mumax)
      } else {
        ## cluster sizes are Poisson conditional on > 0
        csize <- qpois(runif(np, min=dpois(0, mumax)), mumax)
      }
      noff <- sum(csize)
      xparent <- parents$x
      yparent <- parents$y
      x0 <- rep.int(xparent, csize)
      y0 <- rep.int(yparent, csize)
      ## invoke random generator
      dd <- rdisplace(noff, ...)
      mm <- if(is.ppp(dd)) marks(dd) else NULL
      ## validate
      xy <- xy.coords(dd)
      dx <- xy$x
      dy <- xy$y
      if(!(length(dx) == noff))
        stop("rcluster returned the wrong number of points")
      ## create offspring and offspring-to-parent map
      xoff <- x0 + dx
      yoff <- y0 + dy
      parentid <- rep.int(1:np, csize)
      ## trim to window
      retain <- inside.owin(xoff, yoff, win)
      if(is.im(mu))
        retain[retain] <- inside.owin(xoff[retain], yoff[retain], as.owin(mu))
      xoff <- xoff[retain]
      yoff <- yoff[retain]
      parentid <- parentid[retain]
      if(!is.null(mm)) mm <- marksubset(mm, retain)
      ## done
      result <- ppp(xoff, yoff, window=win, check=FALSE, marks=mm)
    }

    if(is.im(mu)) {
      ## inhomogeneously modulated clusters a la Waagepetersen
      P <- eval.im(mu/mumax)
      result <- rthin(result, P)
    }

    if(saveparents) {
      attr(result, "parents") <- parents
      attr(result, "parentid") <- parentid
      attr(result, "expand") <- expand
    }
    
    resultlist[[i]] <- result
  }

  if(nsim == 1 && drop) return(resultlist[[1]])
  names(resultlist) <- paste("Simulation", 1:nsim)
  return(as.solist(resultlist))
}  

rMatClust <- local({
  
  ## like runifdisc but returns only the coordinates
  rundisk <- function(n, radius) {
    R <- radius * sqrt(runif(n, min=0, max=1))
    Theta <- runif(n, min=0, max=2*pi)
    cbind(R * cos(Theta), R * sin(Theta))
  }

  rMatClust <- 
  function(kappa, scale, mu, win = owin(c(0,1),c(0,1)),
           nsim=1, drop=TRUE, saveLambda=FALSE, expand = scale, ...,
           poisthresh=1e-6, saveparents=TRUE) {
    ## Matern Cluster Process with Poisson (mu) offspring distribution
    ## Catch old scale syntax (r)
    if(missing(scale)) scale <- list(...)$r
    check.1.real(scale)
    stopifnot(scale > 0)

    ## trap case of large clusters, close to Poisson
    kok <- is.numeric(kappa) || is.im(kappa)
    if(kok) {
      kappamax <- max(kappa)
    } else {
      kim <- as.im(kappa, W=win, ..., strict=TRUE)
      kra <- range(kim)
      kappamax <- kra[2] + 0.05 * diff(kra)
    }
    if(1/(pi * kappamax * scale^2) < poisthresh) {
      kapmu <- mu * (if(kok) kappa else kim)
      result <- rpoispp(kapmu, win=win, nsim=nsim, drop=drop, warnwin=FALSE)
      return(result)
    }

    result <- rNeymanScott(kappa, scale, list(mu, rundisk), win, radius=scale,
                           nsim=nsim, drop=FALSE,
                           saveparents = saveparents || saveLambda)
    if(saveLambda){
      for(i in 1:nsim) {
        parents <- attr(result[[i]], "parents")
        Lambda <- clusterfield("MatClust", parents, scale=scale, mu=mu, ...)
        attr(result[[i]], "Lambda") <- Lambda[win]
      }
    }
    return(if(nsim == 1 && drop) result[[1]] else result)
  }

  rMatClust
})

                  
rThomas <- local({

  ## random displacements
  gaus <- function(n, sigma) {
    matrix(rnorm(2 * n, mean=0, sd=sigma), ncol=2)
  }

  ## main function
  rThomas <-
      function(kappa, scale, mu, win = owin(c(0,1),c(0,1)), nsim=1, drop=TRUE, 
               saveLambda=FALSE, expand = 4*scale, ...,
               poisthresh=1e-6, saveparents=TRUE) {
      ## Thomas process with Poisson(mu) number of offspring
      ## at isotropic Normal(0,sigma^2) displacements from parent
      ##
      ## Catch old scale syntax (sigma)
      if(missing(scale)) scale <- list(...)$sigma
      check.1.real(scale)
      stopifnot(scale > 0)

      ## trap case of large clusters, close to Poisson
      kok <- is.numeric(kappa) || is.im(kappa)
      if(kok) {
        kappamax <- max(kappa)
      } else {
        kim <- as.im(kappa, W=win, ..., strict=TRUE)
        kra <- range(kim)
        kappamax <- kra[2] + 0.05 * diff(kra)
      }
      if(1/(4*pi * kappamax * scale^2) < poisthresh) {
        kapmu <- mu * (if(kok) kappa else kim)
        result <- rpoispp(kapmu, win=win, nsim=nsim, drop=drop, warnwin=FALSE)
        return(result)
      }
      
      ## determine the maximum radius of clusters
      if(missing(expand))
          expand <- clusterradius("Thomas", scale = scale, ...)

      result <- rNeymanScott(kappa, expand, list(mu, gaus),
                             win, sigma=scale,
                             nsim=nsim, drop=FALSE,
                             saveparents = saveparents || saveLambda)  
      if(saveLambda){
        for(i in 1:nsim) {
          parents <- attr(result[[i]], "parents")
          Lambda <- clusterfield("Thomas", parents, scale=scale, mu=mu, ...)
          attr(result[[i]], "Lambda") <- Lambda[win]
        }
      }
      return(if(nsim == 1 && drop) result[[1]] else result)
    }

  rThomas
})


## ================================================
## Neyman-Scott process with Cauchy kernel function
## ================================================

## scale / omega: scale parameter of Cauchy kernel function
## eta: scale parameter of Cauchy pair correlation function
## eta = 2 * omega

rCauchy <- local({

  ## simulate mixture of normals with inverse-gamma distributed variance
  rnmix.invgam <- function(n = 1, rate) {
    V <- matrix(rnorm(2 * n, 0, 1), nrow = n, ncol = 2)
    s <- 1/rgamma(n, shape=1/2, rate=rate)
    return(sqrt(s) * V)
  }

  ## main function
  rCauchy <- function (kappa, scale, mu, win = owin(), thresh = 0.001,
                       nsim=1, drop=TRUE, saveLambda=FALSE, expand = NULL,
                       ..., poisthresh=1e-6, saveparents=TRUE) {
    ## scale / omega: scale parameter of Cauchy kernel function
    ## eta: scale parameter of Cauchy pair correlation function

    ## Catch old scale syntax (omega)
    dots <- list(...)
    if(missing(scale)) scale <- dots$omega
    
    ## Catch old name 'eps' for 'thresh':
    if(missing(thresh))
        thresh <- dots$eps %orifnull% 0.001

    ## trap case of large clusters, close to Poisson
    kok <- is.numeric(kappa) || is.im(kappa)
    if(kok) {
      kappamax <- max(kappa)
    } else {
      kim <- as.im(kappa, W=win, ..., strict=TRUE)
      kra <- range(kim)
      kappamax <- kra[2] + 0.05 * diff(kra)
    }
    if(1/(pi * kappamax * scale^2) < poisthresh) {
      kapmu <- mu * (if(kok) kappa else kim)
      result <- rpoispp(kapmu, win=win, nsim=nsim, drop=drop, warnwin=FALSE)
      return(result)
    }
    
    ## determine the maximum radius of clusters
    if(missing(expand)){
        expand <- clusterradius("Cauchy", scale = scale, thresh = thresh, ...)
    } else if(!missing(thresh)){
        warning("Argument ", sQuote("thresh"), " is ignored when ", sQuote("expand"), " is given")
    }

    ## simulate
    result <- rNeymanScott(kappa, expand,
                           list(mu, rnmix.invgam),
                           win, rate = scale^2/2, nsim=nsim, drop=FALSE,
                           saveparents = saveparents || saveLambda)
    ## correction from Abdollah: the rate is beta = omega^2 / 2 = eta^2 / 8.
    if(saveLambda){
      for(i in 1:nsim) {
        parents <- attr(result[[i]], "parents")
        Lambda <- clusterfield("Cauchy", parents, scale=scale, mu=mu, ...)
        attr(result[[i]], "Lambda") <- Lambda[win]
      }
    }
    return(if(nsim == 1 && drop) result[[1]] else result)
  }

  rCauchy })

##    
## =================================================================
## Neyman-Scott process with Variance Gamma (Bessel) kernel function
## =================================================================

## nu.ker: smoothness parameter of Variance Gamma kernel function
## omega: scale parameter of kernel function
## nu.pcf: smoothness parameter of Variance Gamma pair correlation function
## eta: scale parameter of Variance Gamma pair correlation function
## nu.pcf = 2 * nu.ker + 1    and    eta = omega

rVarGamma <- local({
  
  ## simulates mixture of isotropic Normal points in 2D with gamma variances
  rnmix.gamma <- function(n = 1, shape, rate) {
    V <- matrix(rnorm(2 * n, 0, 1), nrow = n, ncol = 2)
    s <- rgamma(n, shape=shape, rate=rate)
    return(sqrt(s) * V)
  }

  ## main function
  rVarGamma <- function (kappa, nu, scale, mu, win = owin(),
                         thresh = 0.001, nsim=1, drop=TRUE, saveLambda=FALSE,
                         expand = NULL, ..., poisthresh=1e-6,
                         saveparents=TRUE) {
    ## nu / nu.ker: smoothness parameter of Variance Gamma kernel function
    ## scale / omega: scale parameter of kernel function
    ## Catch old nu.ker/nu.pcf syntax and resolve nu-value.
    dots <- list(...)
    if(missing(nu)){
        nu <- resolve.vargamma.shape(nu.ker=dots$nu.ker, nu.pcf=dots$nu.pcf)$nu.ker
    } else{
        check.1.real(nu)
        stopifnot(nu > -1/2)
    }
    ## Catch old scale syntax (omega)
    if(missing(scale)) scale <- dots$omega
    
    ## Catch old name 'eps' for 'thresh':
    if(missthresh <- missing(thresh))
        thresh <- dots$eps %orifnull% 0.001

    ## trap case of large clusters, close to Poisson
    kok <- is.numeric(kappa) || is.im(kappa)
    if(kok) {
      kappamax <- max(kappa)
    } else {
      kim <- as.im(kappa, W=win, ..., strict=TRUE)
      kra <- range(kim)
      kappamax <- kra[2] + 0.05 * diff(kra)
    }
    if(1/(4 * pi * kappamax * scale^2) < poisthresh) {
      kapmu <- mu * (if(kok) kappa else kim)
      result <- rpoispp(kapmu, win=win, nsim=nsim, drop=drop, warnwin=FALSE)
      return(result)
    }
    
     ## determine the maximum radius of clusters
    if(missing(expand)){
        expand <- clusterradius("VarGamma", scale = scale, nu = nu,
                             thresh = thresh, ...)
    } else if(!missthresh){
        warning("Argument ", sQuote("thresh"), " is ignored when ", sQuote("expand"), " is given")
    }

    ## simulate
    result <- rNeymanScott(kappa, expand,
                           list(mu, rnmix.gamma), win,
##                          WAS:  shape = 2 * (nu.ker + 1)
                           shape = nu + 1,
                           rate = 1/(2 * scale^2),
                           nsim=nsim, drop=FALSE,
                           saveparents = saveparents || saveLambda)
    if(saveLambda){
      for(i in 1:nsim) {
        parents <- attr(result[[i]], "parents")
        Lambda <- clusterfield("VarGamma", parents, scale=scale,
                               nu=nu, mu=mu, ...)
        attr(result[[i]], "Lambda") <- Lambda[win]
      }
    }
    return(if(nsim == 1 && drop) result[[1]] else result)
  }

  rVarGamma
})
