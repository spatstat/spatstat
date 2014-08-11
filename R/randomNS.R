#
#   randomNS.R
#
#   simulating from Neyman-Scott processes
#
#   $Revision: 1.11 $  $Date: 2013/12/10 06:11:43 $
#
#    Original code for rCauchy and rVarGamma by Abdollah Jalilian
#    Other code and modifications by Adrian Baddeley
#    Bug fixes by Abdollah, Adrian, and Rolf Turner

"rNeymanScott" <-
  function(kappa, rmax, rcluster, win = owin(c(0,1),c(0,1)), ..., lmax=NULL)
{
  # Generic Neyman-Scott process
  # Implementation for bounded cluster radius
  #
  # 'rcluster' may be
  #
  #     (1) a function(x,y, ...) that takes the coordinates
  #         (x,y) of the parent point and generates a list(x,y) of offspring
  #
  if(is.function(rcluster))
    return(rPoissonCluster(kappa, rmax, rcluster, win, ..., lmax=lmax))

  #     (2) a list(mu, f) where mu is a numeric value, function, or pixel image
  #         and f is a function(n, ...) generating n i.i.d. offspring at 0,0
  
  if(!(is.list(rcluster) && length(rcluster) == 2))
    stop("rcluster should be either a function, or a list of two elements")
  win <- as.owin(win)
  mu <- rcluster[[1]]
  rdisplace <- rcluster[[2]]
  if(is.numeric(mu)) {
    # homogeneous
    if(!(length(mu) == 1 && mu >= 0))
      stop("rcluster[[1]] should be a single nonnegative number")
    mumax <- mu
  } else if (is.im(mu) || is.function(mu)) {
      # inhomogeneous
    if(is.function(mu)) mu <- as.im(mu, W=win)
    mumax <- max(mu)
  } else stop("rcluster[[1]] should be a number, a function or a pixel image")  
  if(!is.function(rdisplace))
    stop("rcluster[[2]] should be a function")

  # Generate parents in dilated window
  frame <- bounding.box(win)
  dilated <- grow.rectangle(frame, rmax)
  if(is.im(kappa) && !is.subset.owin(dilated, as.owin(kappa)))
    stop(paste("The window in which the image",
               sQuote("kappa"),
               "is defined\n",
               "is not large enough to contain the dilation of the window",
               sQuote("win")))
  parents <- rpoispp(kappa, lmax=lmax, win=dilated)
  np <- npoints(parents)

  # generate cluster sizes
  if(np == 0) {
    # no parents - empty pattern
    result <- ppp(numeric(0), numeric(0), window=win)
    parentid <- integer(0)
  } else {
    csize <- rpois(np, mumax)
    noff <- sum(csize)
    xparent <- parents$x
    yparent <- parents$y
    x0 <- rep.int(xparent, csize)
    y0 <- rep.int(yparent, csize)
    # invoke random generator
    dd <- rdisplace(noff, ...)
    mm <- if(is.ppp(dd)) marks(dd) else NULL
    # validate
    xy <- xy.coords(dd)
    dx <- xy$x
    dy <- xy$y
    if(!(length(dx) == noff))
      stop("rcluster returned the wrong number of points")
    # create offspring and offspring-to-parent map
    xoff <- x0 + dx
    yoff <- y0 + dy
    parentid <- rep.int(1:np, csize)
    # trim to window
    retain <- inside.owin(xoff, yoff, win)
    xoff <- xoff[retain]
    yoff <- yoff[retain]
    parentid <- parentid[retain]
    if(!is.null(mm)) mm <- marksubset(mm, retain)
    # done
    result <- ppp(xoff, yoff, window=win, check=FALSE, marks=mm)
  }

  attr(result, "parents") <- parents
  attr(result, "parentid") <- parentid

  if(is.im(mu)) {
    # inhomogeneously modulated clusters a la Waagepetersen
    P <- eval.im(mu/mumax)
    result <- rthin(result, P)
  }
  return(result)
}  

"rMatClust" <- local({
  
  # like runifdisc but returns only the coordinates
  rundisk <- function(n, radius) {
    R <- radius * sqrt(runif(n, min=0, max=1))
    Theta <- runif(n, min=0, max=2*pi)
    cbind(R * cos(Theta), R * sin(Theta))
  }

  rMatClust <- 
  function(kappa, r, mu, win = owin(c(0,1),c(0,1))) {
    # Matern Cluster Process with Poisson (mu) offspring distribution
    stopifnot(is.numeric(r) && length(r) == 1 && r > 0)
    result <- rNeymanScott(kappa, r, list(mu, rundisk), win, radius=r)  
    return(result)
  }

  rMatClust
})

                  
"rThomas" <- local({

  # random displacements
  gaus <- function(n, sigma) {
    matrix(rnorm(2 * n, mean=0, sd=sigma), ncol=2)
  }

  # main function
  rThomas <-
    function(kappa, sigma, mu, win = owin(c(0,1),c(0,1))) {
      # Thomas process with Poisson(mu) number of offspring
      # at isotropic Normal(0,sigma^2) displacements from parent
      #
      stopifnot(is.numeric(sigma) && length(sigma) == 1 && sigma > 0)

      result <- rNeymanScott(kappa, 4 * sigma, list(mu, gaus),
                             win, sigma=sigma)  
      return(result)
    }
  rThomas
})


# ================================================
# Neyman-Scott process with Cauchy kernel function
# ================================================

# omega: scale parameter of Cauchy kernel function
# eta: scale parameter of Cauchy pair correlation function
# eta = 2 * omega

rCauchy <- local({

  # simulate mixture of normals with inverse-gamma distributed variance
  rnmix.invgam <- function(n = 1, rate) {
    V <- matrix(rnorm(2 * n, 0, 1), nrow = n, ncol = 2)
    s <- 1/rgamma(n, shape=1/2, rate=rate)
    return(sqrt(s) * V)
  }

  # threshold the kernel function in polar coordinate
  kernthresh <- function(r, eta, eps) {
    4 * (r/eta^2)/((1 + (2 * r/eta)^2)^(3/2)) - eps
  }
  
  # main function
  rCauchy <- function (kappa, omega, mu, win = owin(), eps = 0.001) {

    # omega: scale parameter of Cauchy kernel function
    # eta: scale parameter of Cauchy pair correlation function
    eta     <- 2 * omega
    
    # determine the maximum radius of clusters
    rmax <- uniroot(kernthresh,
                    lower = eta/2, upper = 5 * diameter(as.rectangle(win)),
                    eta = eta, eps = eps)$root
    # simulate
    result <- rNeymanScott(kappa, rmax,
                           list(mu, rnmix.invgam),
                           win, rate = eta^2/8)
    # correction from Abdollah: the rate is beta = omega^2 / 2 = eta^2 / 8.
    return(result)
  }

  rCauchy })

#    
# =================================================================
# Neyman-Scott process with Variance Gamma (Bessel) kernel function
# =================================================================

# nu.ker: smoothness parameter of Variance Gamma kernel function
# omega: scale parameter of kernel function
# nu.pcf: smoothness parameter of Variance Gamma pair correlation function
# eta: scale parameter of Variance Gamma pair correlation function
# nu.pcf = 2 * nu.ker + 1    and    eta = omega

rVarGamma <- local({
  
  # simulates mixture of isotropic Normal points in 2D with gamma variances
  rnmix.gamma <- function(n = 1, shape, rate) {
    V <- matrix(rnorm(2 * n, 0, 1), nrow = n, ncol = 2)
    s <- rgamma(n, shape=shape, rate=rate)
    return(sqrt(s) * V)
  }

  # kernel function in polar coordinates
  kernfun.old <- function(r, nu.ker, omega, eps) {
    numer <- ((r/omega)^(nu.ker+1)) * besselK(r/omega, nu.ker)
    denom <- (2^nu.ker) * omega * gamma(nu.ker + 1)
    numer/denom - eps
  }
  kernfun <- function(r, nu.ker, omega, eps) {
    numer <- ((r/omega)^(nu.ker + 1)) * besselK(r/omega, nu.ker)
    denom <- pi * (2^(nu.ker+1)) * omega^2 * gamma(nu.ker + 1)
    numer/denom - eps
  }
  
  # main function
  rVarGamma <- function (kappa, nu.ker=NULL, omega, mu, win = owin(),
                         eps = 0.001, nu.pcf=NULL) {
    # nu.ker: smoothness parameter of Variance Gamma kernel function
    # omega: scale parameter of kernel function

    nu.ker <- resolve.vargamma.shape(nu.ker=nu.ker, nu.pcf=nu.pcf)$nu.ker
    
    # determine the maximum radius of clusters
    rmax <- uniroot(kernfun,
                    lower = omega, upper = 5 * diameter(as.rectangle(win)),
                    nu.ker = nu.ker, omega=omega, eps=eps)$root
    # simulate
    result <- rNeymanScott(kappa, rmax,
                           list(mu, rnmix.gamma), win,
#                          WAS:  shape = 2 * (nu.ker + 1)
                           shape = nu.ker + 1,
                           rate = 1/(2 * omega^2))
    return(result)
  }

  rVarGamma })

