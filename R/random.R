##
##    random.R
##
##    Functions for generating random point patterns
##
##    $Revision: 4.76 $   $Date: 2015/03/13 11:15:51 $
##
##
##    runifpoint()      n i.i.d. uniform random points ("binomial process")
##
##    runifpoispp()     uniform Poisson point process
##
##    rpoispp()         general Poisson point process (thinning method)
##
##    rpoint()          n independent random points (rejection/pixel list)
##
##    rMaternI()        Mat'ern model I 
##    rMaternII()       Mat'ern model II
##    rSSI()            Simple Sequential Inhibition process
##
##    rthin()           independent random thinning
##    rjitter()         random perturbation
##
##    Examples:
##          u01 <- owin(0:1,0:1)
##          plot(runifpoispp(100, u01))
##          X <- rpoispp(function(x,y) {100 * (1-x/2)}, 100, u01)
##          X <- rpoispp(function(x,y) {ifelse(x < 0.5, 100, 20)}, 100)
##          plot(X)
##          plot(rMaternI(100, 0.02))
##          plot(rMaternII(100, 0.05))
##

runifrect <- function(n, win=owin(c(0,1),c(0,1)), nsim=1, drop=TRUE)
{
  ## no checking
  xr <- win$xrange
  yr <- win$yrange
  if(nsim == 1) {
    x <- runif(n, min=xr[1], max=xr[2])
    y <- runif(n, min=yr[1], max=yr[2])
    X <- ppp(x, y, window=win, check=FALSE)
    return(if(drop) X else solist(X))
  } 
  result <- vector(mode="list", length=nsim)
  for(isim in 1:nsim) {
    x <- runif(n, min=xr[1], max=xr[2])
    y <- runif(n, min=yr[1], max=yr[2])
    result[[isim]] <- ppp(x, y, window=win, check=FALSE)
  }
  return(as.solist(result))
}

runifdisc <- function(n, radius=1, centre=c(0,0), ..., nsim=1, drop=TRUE)
{
  ## i.i.d. uniform points in the disc of radius r and centre (x,y)
  disque <- disc(centre=centre, radius=radius, ...)
  if(nsim == 1) {
    theta <- runif(n, min=0, max= 2 * pi)
    s <- sqrt(runif(n, min=0, max=radius^2))
    X <- ppp(centre[1] + s * cos(theta),
             centre[2] + s * sin(theta),
             window=disque, check=FALSE)
    return(if(drop) X else solist(X))
  }
  result <- vector(mode="list", length=nsim)
  for(isim in 1:nsim) {
    theta <- runif(n, min=0, max= 2 * pi)
    s <- sqrt(runif(n, min=0, max=radius^2))
    result[[isim]] <- ppp(centre[1] + s * cos(theta),
                          centre[2] + s * sin(theta),
                          window=disque, check=FALSE)
  }
  return(as.solist(result))
}


runifpoint <- function(n, win=owin(c(0,1),c(0,1)),
                       giveup=1000, warn=TRUE, ...,
                       nsim=1, drop=TRUE, ex=NULL)
{
  if(missing(n) && missing(win) && !is.null(ex)) {
    stopifnot(is.ppp(ex))
    n <- npoints(ex)
    win <- Window(ex)
  } else {
    win <- as.owin(win)
    check.1.integer(n)
    stopifnot(n >= 0)
  }

  if(n == 0) {
    emp <- ppp(numeric(0), numeric(0), window=win)
    if(nsim == 1) return(emp)
    result <- rep(list(emp), nsim)
    names(result) <- paste("Simulation", 1:nsim)
    return(as.solist(result))
  }

  if(warn) {
    nhuge <- spatstat.options("huge.npoints")
    if(n > nhuge) {
      whinge <- paste("Attempting to generate", n, "random points")
      message(whinge)
      warning(whinge, call.=FALSE)
    }
  }

  switch(win$type,
         rectangle = {
           return(runifrect(n, win, nsim=nsim))
         },
         mask = {
           dx <- win$xstep
           dy <- win$ystep
           ## extract pixel coordinates and probabilities
           rxy <- rasterxy.mask(win, drop=TRUE)
           xpix <- rxy$x
           ypix <- rxy$y
           ## make a list of nsim point patterns
           result <- vector(mode="list", length=nsim)
           for(isim in 1:nsim) {
             ## select pixels with equal probability
             id <- sample(seq_along(xpix), n, replace=TRUE)
             ## extract pixel centres and randomise within pixels
             x <- xpix[id] + runif(n, min= -dx/2, max=dx/2)
             y <- ypix[id] + runif(n, min= -dy/2, max=dy/2)
             result[[isim]] <- ppp(x, y, window=win, check=FALSE)
           }
         },
         polygonal={
           ## make a list of nsim point patterns
           result <- vector(mode="list", length=nsim)
           for(isim in 1:nsim) {
             ## rejection method
             ## initialise empty pattern
             x <- numeric(0)
             y <- numeric(0)
             X <- ppp(x, y, window=win)
             ##
             ## rectangle in which trial points will be generated
             box <- boundingbox(win)
             ## 
             ntries <- 0
             repeat {
               ntries <- ntries + 1
               ## generate trial points in batches of n
               qq <- runifrect(n, box) 
               ## retain those which are inside 'win'
               qq <- qq[win]
               ## add them to result
               X <- superimpose(X, qq, W=win, check=FALSE)
               ## if we have enough points, exit
               if(X$n > n) {
                 result[[isim]] <- X[1:n]
                 break
               } else if(X$n == n) {
                 result[[isim]] <- X
                 break
               } else if(ntries >= giveup) {
                 ## otherwise get bored eventually
                 stop(paste("Gave up after", giveup * n, "trials,",
                            X$n, "points accepted"))
               }
             }
           }
         },
         stop("Unrecognised window type")
         )
  
  ## list of point patterns produced.
  if(nsim == 1 && drop)
    return(result[[1]])
  names(result) <- paste("Simulation", 1:nsim)
  return(as.solist(result))
}

runifpoispp <- function(lambda, win = owin(c(0,1),c(0,1)), ...,
                        nsim=1, drop=TRUE) {
  win <- as.owin(win)
  if(!is.numeric(lambda) || length(lambda) > 1 ||
     !is.finite(lambda) || lambda < 0)
    stop("Intensity lambda must be a single finite number >= 0")

  if(lambda == 0) {
    ## return empty pattern
    emp <- ppp(numeric(0), numeric(0), window=win)
    if(nsim == 1 && drop) return(emp)
    result <- rep(list(emp), nsim)
    names(result) <- paste("Simulation", 1:nsim)
    return(as.solist(result))
  }

  ## will generate Poisson process in enclosing rectangle and trim it
  box <- boundingbox(win)
  meanN <- lambda * area(box)
  
  if(nsim == 1) {
    n <- rpois(1, meanN)
    if(!is.finite(n))
      stop(paste("Unable to generate Poisson process with a mean of",
                 meanN, "points"))
    X <- runifpoint(n, box)
    ## trim to window
    if(win$type != "rectangle")
      X <- X[win]
    return(if(drop) X else solist(X))
  }
  result <- vector(mode="list", length=nsim)
  for(isim in 1:nsim) {
    n <- rpois(1, meanN)
    if(!is.finite(n))
      stop(paste("Unable to generate Poisson process with a mean of",
                 meanN, "points"))
    X <- runifpoint(n, box)
    ## trim to window
    if(win$type != "rectangle")
      X <- X[win]
    result[[isim]] <- X
  }
  names(result) <- paste("Simulation", 1:nsim)
  return(as.solist(result))
}

rpoint <- function(n, f, fmax=NULL,
                   win=unit.square(), ..., giveup=1000,verbose=FALSE,
                   nsim=1, drop=TRUE) {
  
  if(missing(f) || (is.numeric(f) && length(f) == 1))
    ## uniform distribution
    return(runifpoint(n, win, giveup, nsim=nsim, drop=drop))
  
  ## non-uniform distribution....
  if(!is.function(f) && !is.im(f))
    stop(paste(sQuote("f"),
               "must be either a function or an",
               sQuote("im"), "object"))
  
  if(is.im(f)) {
    ## ------------ PIXEL IMAGE ---------------------
    wf <- as.owin(f)
    if(n == 0) {
      ## return empty pattern(s)
      emp <- ppp(numeric(0), numeric(0), window=wf)
      if(nsim == 1 && drop) return(emp)
      result <- rep(list(emp), nsim)
      names(result) <- paste("Simulation", 1:nsim)
      return(as.solist(result))
    }
    w <- as.mask(wf)
    M <- w$m
    dx <- w$xstep
    dy <- w$ystep
    ## extract pixel coordinates and probabilities
    rxy <- rasterxy.mask(w, drop=TRUE)
    xpix <- rxy$x
    ypix <- rxy$y
    ppix <- as.vector(f$v[M]) ## not normalised - OK
    ##
    if(nsim == 1) {
      ## select pixels
      id <- sample(length(xpix), n, replace=TRUE, prob=ppix)
      ## extract pixel centres and randomise within pixels
      x <- xpix[id] + runif(n, min= -dx/2, max=dx/2)
      y <- ypix[id] + runif(n, min= -dy/2, max=dy/2)
      X <- ppp(x, y, window=wf, check=FALSE)
      return(if(drop) X else solist(X))
    }
    result <- vector(mode="list", length=nsim)
    for(isim in 1:nsim) {
      ## select pixels
      id <- sample(length(xpix), n, replace=TRUE, prob=ppix)
      ## extract pixel centres and randomise within pixels
      x <- xpix[id] + runif(n, min= -dx/2, max=dx/2)
      y <- ypix[id] + runif(n, min= -dy/2, max=dy/2)
      result[[isim]] <- ppp(x, y, window=wf, check=FALSE)
    }
    names(result) <- paste("Simulation", 1:nsim)
    return(as.solist(result))
  }

  ## ------------ FUNCTION  ---------------------  
  ## Establish parameters for rejection method

  verifyclass(win, "owin")

  if(n == 0) {
    ## return empty pattern(s)
    emp <- ppp(numeric(0), numeric(0), window=win)
    if(nsim == 1 && drop) return(emp)
    result <- rep(list(emp), nsim)
    names(result) <- paste("Simulation", 1:nsim)
    return(as.solist(result))
  }
  
  if(is.null(fmax)) {
    ## compute approx maximum value of f
    imag <- as.im(f, win, ...)
    summ <- summary(imag)
    fmax <- summ$max + 0.05 * diff(summ$range)
  }
  irregular <- (win$type != "rectangle")
  box <- boundingbox(win)

  result <- vector(mode="list", length=nsim)
  for(isim in 1:nsim) {

    ## initialise empty pattern
    X <- ppp(numeric(0), numeric(0), window=win)
  
    pbar <- 1
    nremaining <- n
    totngen <- 0
    
    ## generate uniform random points in batches
    ## and apply the rejection method.
    ## Collect any points that are retained in X

    ntries <- 0
    repeat{
      ntries <- ntries + 1
      ## proposal points
      ngen <- nremaining/pbar + 10
      totngen <- totngen + ngen
      prop <- runifrect(ngen, box)
      if(irregular)
        prop <- prop[win]
      if(prop$n > 0) {
        fvalues <- f(prop$x, prop$y, ...)
        paccept <- fvalues/fmax
        u <- runif(prop$n)
        ## accepted points
        Y <- prop[u < paccept]
        if(Y$n > 0) {
          ## add to X
          X <- superimpose(X, Y, W=win, check=FALSE)
          nX <- X$n
          pbar <- nX/totngen
          nremaining <- n - nX
          if(nremaining <= 0) {
            ## we have enough!
            if(verbose)
              splat("acceptance rate = ", round(100 * pbar, 2), "%")
            result[[isim]] <- if(nX == n) X else X[1:n]
            break
          }
        }
      }
      if(ntries > giveup)
        stop(paste("Gave up after",giveup * n,"trials with",
                   X$n, "points accepted"))
    }
  }
  if(nsim == 1 && drop) return(result[[1]])
  names(result) <- paste("Simulation", 1:nsim)
  return(as.solist(result))
}

rpoispp <- function(lambda, lmax=NULL, win = owin(), ...,
                    nsim=1, drop=TRUE, ex=NULL) {
  ## arguments:
  ##     lambda  intensity: constant, function(x,y,...) or image
  ##     lmax     maximum possible value of lambda(x,y,...)
  ##     win     default observation window (of class 'owin')
  ##   ...       arguments passed to lambda(x, y, ...)
  ##     nsim    number of replicate simulations

  if(missing(lambda) && is.null(lmax) && missing(win) && !is.null(ex)) {
    lambda <- intensity(unmark(ex))
    win <- Window(ex)
  } else {
    if(!(is.numeric(lambda) || is.function(lambda) || is.im(lambda)))
      stop(paste(sQuote("lambda"),
                 "must be a constant, a function or an image"))
    if(is.numeric(lambda) && !(length(lambda) == 1 && lambda >= 0))
      stop(paste(sQuote("lambda"),
                 "must be a single, nonnegative number"))
    if(!is.null(lmax)) {
      if(!is.numeric(lmax))
        stop("lmax should be a number")
      if(length(lmax) > 1)
        stop("lmax should be a single number")
    }
    win <- if(is.im(lambda))
      rescue.rectangle(as.owin(lambda))
    else
      as.owin(win)
  }
  
  if(is.numeric(lambda)) 
    ## uniform Poisson
    return(runifpoispp(lambda, win, nsim=nsim, drop=drop))

  ## inhomogeneous Poisson
  ## perform thinning of uniform Poisson

  if(is.null(lmax)) {
    imag <- as.im(lambda, win, ...)
    summ <- summary(imag)
    lmax <- summ$max + 0.05 * diff(summ$range)
  } 

  if(is.function(lambda)) {
    ## function lambda
    if(nsim == 1) {
      X <- runifpoispp(lmax, win)  ## includes sanity checks on `lmax'
      if(X$n == 0) return(X)
      prob <- lambda(X$x, X$y, ...)/lmax
      u <- runif(X$n)
      retain <- (u <= prob)
      X <- X[retain, ]
      return(if(drop) X else solist(X))
    }
    result <- runifpoispp(lmax, win, nsim=nsim, drop=FALSE)
    for(isim in 1:nsim) {
      X <- result[[isim]]
      if(X$n > 0) {
        prob <- lambda(X$x, X$y, ...)/lmax
        u <- runif(X$n)
        retain <- (u <= prob)
        result[[isim]] <- X[retain, ]
      }
    }
    return(as.solist(result))
  }

  if(is.im(lambda)) {
    ## image lambda
    if(nsim == 1) {
      X <- runifpoispp(lmax, win)
      if(X$n == 0) return(X)
      prob <- lambda[X]/lmax
      u <- runif(X$n)
      retain <- (u <= prob)
      X <- X[retain, ]
      return(if(drop) X else solist(X))
    }
    result <- runifpoispp(lmax, win, nsim=nsim, drop=FALSE)
    for(isim in 1:nsim) {
      X <- result[[isim]]
      if(X$n > 0) {
        prob <- lambda[X]/lmax
        u <- runif(X$n)
        retain <- (u <= prob)
        result[[isim]] <- X[retain, ]
      }
    }
    return(as.solist(result))
  }
  stop(paste(sQuote("lambda"), "must be a constant, a function or an image"))
}
    
rMaternI <- function(kappa, r, win = owin(c(0,1),c(0,1)), stationary=TRUE,
                     ..., nsim=1, drop=TRUE)
{
  win <- as.owin(win)
  stopifnot(is.numeric(r) && length(r) == 1)
  if(stationary) {
    ## generate in a larger window
    bigbox <- grow.rectangle(as.rectangle(win), r)
    X <- rpoispp(kappa, win=bigbox, nsim=nsim)
  } else {
    X <- rpoispp(kappa, win=win, nsim=nsim)
  }
  result <- if(nsim == 1) solist(X) else X
  for(isim in 1:nsim) {
    Y <- result[[isim]]
    if(npoints(Y) > 1) {
      d <- nndist(Y)
      Y <- Y[d > r]
    }
    if(stationary)
      Y <- Y[win]
    result[[isim]] <- Y
  }
  if(nsim == 1 && drop) return(result[[1]])
  return(as.solist(result))
}
    
rMaternII <- function(kappa, r, win = owin(c(0,1),c(0,1)), stationary=TRUE,
                      ..., nsim=1, drop=TRUE)
{
  win <- as.owin(win)
  stopifnot(is.numeric(r) && length(r) == 1)
  if(stationary) {
    bigbox <- grow.rectangle(as.rectangle(win), r)
    X <- rpoispp(kappa, win=bigbox, nsim=nsim)
  } else {
    X <- rpoispp(kappa, win=win, nsim=nsim)
  }
  result <- if(nsim == 1) list(X) else X
  for(isim in 1:nsim) {
    Y <- result[[isim]]
    nY <- npoints(Y)
    if(nY > 1) {
      ## matrix of squared pairwise distances
      d2 <- pairdist(Y, squared=TRUE)
      close <- (d2 <= r^2)
      ## random order 1:n
      age <- sample(seq_len(nY), nY, replace=FALSE)
      earlier <- outer(age, age, ">")
      conflict <- close & earlier
      ## delete <- apply(conflict, 1, any)
      delete <- matrowany(conflict)
      Y <- Y[!delete]
    }
    if(stationary)
      Y <- Y[win]
    result[[isim]] <- Y
  }
  if(nsim == 1 && drop) return(result[[1]])
  return(as.solist(result))
}
  
rSSI <- function(r, n=Inf, win = square(1), 
                 giveup = 1000, x.init=NULL, ...,
                 f=NULL, fmax=NULL,
                 nsim=1, drop=TRUE)
{
  win.given <- !missing(win) && !is.null(win)
  stopifnot(is.numeric(r) && length(r) == 1 && r >= 0)
  stopifnot(is.numeric(n) && length(n) == 1 && n >= 0)
  ##
  if(!is.null(f)) {
    stopifnot(is.numeric(f) || is.im(f) || is.function(f))
    if(is.null(fmax) && !is.numeric(f))
      fmax <- if(is.im(f)) max(f) else max(as.im(f, win))
  }
  ##
  if(nsim > 1) {
    result <- vector(mode="list", length=nsim)
    if(!win.given) win <- square(1)
    for(isim in 1:nsim) {
      progressreport(isim, nsim)
      result[[isim]] <- rSSI(r=r, n=n, win=win, giveup=giveup, x.init=x.init,
                             f=f, fmax=fmax)
    }
    names(result) <- paste("Simulation", 1:nsim)
    return(as.solist(result))
  }
  ## Simple Sequential Inhibition process
  ## fixed number of points
  ## Naive implementation, proposals are uniform
  if(is.null(x.init)) {
    ## start with empty pattern in specified window
    win <- as.owin(win)
    x.init <- ppp(numeric(0),numeric(0), window=win)
  } else {
    ## start with specified pattern
    stopifnot(is.ppp(x.init))
    if(!win.given) {
      win <- as.owin(x.init)
    } else {
      ## check compatibility of windows
      if(!identical(win, as.owin(x.init)))
        warning(paste("Argument", sQuote("win"),
                      "is not the same as the window of", sQuote("x.init")))
      x.init.new <- x.init[win]
      if(npoints(x.init.new) == 0)
        stop(paste("No points of x.init lie inside the specified window",
                   sQuote("win")))
      nlost <- npoints(x.init) - npoints(x.init.new)
      if(nlost > 0) 
        warning(paste(nlost, "out of",
                      npoints(x.init), "points of the pattern x.init",
                      "lay outside the specified window",
                      sQuote("win")))
      x.init <- x.init.new
    }
    if(n < npoints(x.init))
      stop(paste("x.init contains", npoints(x.init), "points",
                 "but a pattern containing only n =", n, "points", 
                 "is required"))
    if(n == npoints(x.init)) {
      warning(paste("Initial state x.init already contains", n, "points;",
                    "no further points were added"))
      return(if(drop) x.init else solist(x.init))
    }
  }
  X <- x.init
  r2 <- r^2
  if(!is.infinite(n) && (n * pi * r2/4  > area(win)))
    warning(paste("Window is too small to fit", n, "points",
               "at minimum separation", r))
  ntries <- 0
  while(ntries < giveup) {
    ntries <- ntries + 1
    qq <- if(is.null(f)) runifpoint(1, win) else rpoint(1, f, fmax, win)
    dx <- qq$x[1] - X$x
    dy <- qq$y[1] - X$y
    if(all(dx^2 + dy^2 > r2))
      X <- superimpose(X, qq, W=win, check=FALSE)
    if(X$n >= n)
      return(if(drop) X else solist(X))
  }
  if(!is.infinite(n))
    warning(paste("Gave up after", giveup,
                  "attempts with only", X$n, "points placed out of", n))
  return(if(drop) X else solist(X))
}

rPoissonCluster <-
  function(kappa, expand, rcluster, win = owin(c(0,1),c(0,1)), ...,
           lmax=NULL, nsim=1, drop=TRUE)
{
  ## Generic Poisson cluster process
  ## Implementation for bounded cluster radius
  ##
  ## 'rcluster' is a function(x,y) that takes the coordinates
  ## (x,y) of the parent point and generates a list(x,y) of offspring
  ##
  ## "..." are arguments to be passed to 'rcluster()'
  ##

  ## Catch old argument name rmax for expand, and allow rmax to be
  ## passed to rcluster (and then be ignored)
  if(missing(expand) && !is.null(rmax <- list(...)$rmax)){
      expand <- rmax
      f <- rcluster
      rcluster <- function(..., rmax) f(...)
  }
  win <- as.owin(win)
  
  ## Generate parents in dilated window
  frame <- boundingbox(win)
  dilated <- owin(frame$xrange + c(-expand, expand),
                  frame$yrange + c(-expand, expand))
  if(is.im(kappa) && !is.subset.owin(dilated, as.owin(kappa)))
    stop(paste("The window in which the image",
               sQuote("kappa"),
               "is defined\n",
               "is not large enough to contain the dilation of the window",
               sQuote("win")))
  parentlist <- rpoispp(kappa, lmax=lmax, win=dilated, nsim=nsim)
  if(nsim == 1) parentlist <- list(parentlist)

  resultlist <- vector(mode="list", length=nsim)
  for(isim in 1:nsim) {
    parents <- parentlist[[isim]]
    result <- NULL
    ## generate clusters
    np <- parents$n
    if(np > 0) {
      xparent <- parents$x
      yparent <- parents$y
      for(i in seq_len(np)) {
        ## generate random offspring of i-th parent point
        cluster <- rcluster(xparent[i], yparent[i], ...)
        if(!inherits(cluster, "ppp"))
          cluster <- ppp(cluster$x, cluster$y, window=frame, check=FALSE)
        ## skip if cluster is empty
        if(cluster$n > 0) {
          ## trim to window
          cluster <- cluster[win]
          if(is.null(result)) {
            ## initialise offspring pattern and offspring-to-parent map
            result <- cluster
            parentid <- rep.int(1, cluster$n)
          } else {
            ## add to pattern
            result <- superimpose(result, cluster, W=win, check=FALSE)
            ## update offspring-to-parent map
            parentid <- c(parentid, rep.int(i, cluster$n))
          }
        }
      }
    } else {
      ## no parents - empty pattern
      result <- ppp(numeric(0), numeric(0), window=win)
      parentid <- integer(0)
    }

    attr(result, "parents") <- parents
    attr(result, "parentid") <- parentid
    attr(result, "expand") <- expand

    resultlist[[isim]] <- result
  }

  if(nsim == 1 && drop) return(resultlist[[1]])

  names(resultlist) <- paste("Simulation", 1:nsim)
  return(as.solist(resultlist))
}  

rGaussPoisson <- local({
  
  rGaussPoisson <- function(kappa, r, p2, win=owin(c(0,1), c(0,1)),
                            ..., nsim=1, drop=TRUE) {
    ## Gauss-Poisson process
    result <- rPoissonCluster(kappa, 1.05 * r, oneortwo,
                              win, radius=r/2, p2=p2, nsim=nsim, drop=drop)
    return(result)
  }

  oneortwo <- function(x0, y0, radius, p2) {
    if(runif(1) > p2) 
      ## one point
      return(list(x=x0, y=y0))
    ## two points
    theta <- runif(1, min=0, max=2*pi)
    return(list(x=x0+c(-1,1)*radius*cos(theta),
                y=y0+c(-1,1)*radius*sin(theta)))
  }

  rGaussPoisson
})

  
rstrat <- function(win=square(1), nx, ny=nx, k=1, nsim=1, drop=TRUE) {
  win <- as.owin(win)
  stopifnot(nx >= 1 && ny >= 1)
  stopifnot(k >= 1)
  if(nsim == 1) {
    xy <- stratrand(win, nx, ny, k)
    Xbox <- ppp(xy$x, xy$y, win$xrange, win$yrange, check=FALSE)
    X <- Xbox[win]
    return(if(drop) X else solist(X))
  }
  result <- vector(mode="list", length=nsim)
  for(isim in 1:nsim) {
    xy <- stratrand(win, nx, ny, k)
    Xbox <- ppp(xy$x, xy$y, win$xrange, win$yrange, check=FALSE)
    result[[isim]] <- Xbox[win]
  }
  names(result) <- paste("Simulation", 1:nsim)
  return(as.solist(result))
}

xy.grid <- function(xr, yr, nx, ny, dx, dy) {
  nx.given <- !is.null(nx)
  ny.given <- !is.null(ny)
  dx.given <- !is.null(dx)
  dy.given <- !is.null(dy)
  if(nx.given && dx.given)
    stop("Do not give both nx and dx")    
  if(nx.given) {
    stopifnot(nx >= 1)
    x0 <- seq(from=xr[1], to=xr[2], length.out=nx+1)
    dx <- diff(xr)/nx
  } else if(dx.given) {
    stopifnot(dx > 0)
    x0 <- seq(from=xr[1], to=xr[2], by=dx)
    nx <- length(x0) - 1
  } else stop("Need either nx or dx")
  ## determine y grid
  if(ny.given && dy.given)
    stop("Do not give both ny and dy")    
  if(ny.given) {
    stopifnot(ny >= 1)
    y0 <- seq(from=yr[1], to=yr[2], length.out=ny+1)
    dy <- diff(yr)/ny
  } else {
    if(is.null(dy)) dy <- dx
    stopifnot(dy > 0)
    y0 <- seq(from=yr[1], to=yr[2], by=dy)
    ny <- length(y0) - 1
  }
  return(list(x0=x0, y0=y0, nx=nx, ny=ny, dx=dx, dy=dy))
}
  
rsyst <- function(win=square(1), nx=NULL, ny=nx, ..., dx=NULL, dy=dx,
                  nsim=1, drop=TRUE) {
  win <- as.owin(win)
  xr <- win$xrange
  yr <- win$yrange
  ## determine grid coordinates 
  if(missing(ny)) ny <- NULL
  if(missing(dy)) dy <- NULL
  g <- xy.grid(xr, yr, nx, ny, dx, dy)
  x0 <- g$x0
  y0 <- g$y0
  dx <- g$dx
  dy <- g$dy
  ## assemble grid and randomise location
  xy0 <- expand.grid(x=x0, y=y0)
  if(nsim == 1) {
    x <- xy0$x + runif(1, min = 0, max = dx)
    y <- xy0$y + runif(1, min = 0, max = dy)
    Xbox <- ppp(x, y, xr, yr, check=FALSE)
    ## trim to window
    X <- Xbox[win]
    return(if(drop) X else solist(X))
  }
  result <- vector(mode="list", length=nsim)
  for(isim in 1:nsim) {
    x <- xy0$x + runif(1, min = 0, max = dx)
    y <- xy0$y + runif(1, min = 0, max = dy)
    Xbox <- ppp(x, y, xr, yr, check=FALSE)
    ## trim to window
    result[[isim]] <- Xbox[win]
  }
  names(result) <- paste("Simulation", 1:nsim)
  return(as.solist(result))
}

rcellnumber <- function(n, N=10) {
  if(!missing(N)) {
    if(round(N) != N) stop("N must be an integer")
    stopifnot(is.finite(N))
    stopifnot(N > 1)
  }
  u <- runif(n, min=0, max=1)
  p0 <- 1/N
  pN <- 1/(N * (N-1))
  k <- ifelse(u < p0, 0, ifelse(u < (1 - pN), 1, N))
  return(k)
}

rcell <- function(win=square(1), nx=NULL, ny=nx, ...,
                  dx=NULL, dy=dx, N=10, nsim=1, drop=TRUE) {
  win <- as.owin(win)
  xr <- win$xrange
  yr <- win$yrange
  ## determine grid coordinates 
  if(missing(ny)) ny <- NULL
  if(missing(dy)) dy <- NULL
  g <- xy.grid(xr, yr, nx, ny, dx, dy)
  nx <- g$nx
  ny <- g$ny
  x0 <- g$x0
  y0 <- g$y0
  dx <- g$dx
  dy <- g$dy
  ## generate pattern(s)
  result <- vector(mode="list", length=nsim)
  for(isim in 1:nsim) {
    x <- numeric(0)
    y <- numeric(0)
    for(ix in seq_len(nx))
      for(iy in seq_len(ny)) {
        nij <- rcellnumber(1, N)
        x <- c(x, x0[ix] + runif(nij, min=0, max=dx))
        y <- c(y, y0[iy] + runif(nij, min=0, max=dy))
      }
    Xbox <- ppp(x, y, xr, yr, check=FALSE)
    result[[isim]] <- Xbox[win]
  }
  if(nsim == 1 && drop) return(result[[1]])
  names(result) <- paste("Simulation", 1:nsim)
  return(as.solist(result))
}


rthin <- function(X, P, ..., nsim=1, drop=TRUE) {
  verifyclass(X, "ppp")

  nX <- npoints(X)
  if(nX == 0) {
    if(nsim == 1 && drop) return(X)
    result <- rep(list(X), nsim)
    names(result) <- paste("Simulation", 1:nsim)
    return(as.solist(result))
  }

  if(is.numeric(P)) {
    ## vector of retention probabilities
    pX <- P
    if(length(pX) != nX) {
      if(length(pX) == 1)
        pX <- rep.int(pX, nX)
      else 
        stop("Length of vector P does not match number of points of X")
    }
    if(any(is.na(pX)))
      stop("P contains NA's")
  } else if(is.function(P)) {
    ## function - evaluate it at points of X
    pX <- P(X$x, X$y, ...)
    if(length(pX) != nX)
      stop("Function P returned a vector of incorrect length")
    if(!is.numeric(pX))
      stop("Function P returned non-numeric values")
    if(any(is.na(pX)))
      stop("Function P returned some NA values")
  } else if(is.im(P)) {
    ## image - look it up
    if(!(P$type %in% c("integer", "real")))
      stop("Values of image P should be numeric")
    pX <- P[X, drop=FALSE]
    if(any(is.na(pX)))
      stop("some points of X lie outside the domain of image P")
  } else
  stop("Unrecognised format for P")

  if(min(pX) < 0) stop("some probabilities are negative")
  if(max(pX) > 1) stop("some probabilities are greater than 1")

  if(nsim == 1) {
    retain <- (runif(length(pX)) < pX)
    Y <- X[retain]
    ## also handle offspring-to-parent map if present
    if(!is.null(parentid <- attr(X, "parentid")))
      attr(Y, "parentid") <- parentid[retain]
    return(if(drop) Y else solist(Y))
  }
  result <- vector(mode="list", length=nsim)
  for(isim in 1:nsim) {
    retain <- (runif(length(pX)) < pX)
    Y <- X[retain]
    ## also handle offspring-to-parent map if present
    if(!is.null(parentid <- attr(X, "parentid")))
      attr(Y, "parentid") <- parentid[retain]
    result[[isim]] <- Y
  }
  names(result) <- paste("Simulation", 1:nsim)
  return(as.solist(result))
}


## rjitter

rjitter <- function(X, radius, retry=TRUE, giveup=10000, ...,
                    nsim=1, drop=TRUE) {
  verifyclass(X, "ppp")
  if(missing(radius) || is.null(radius))
    radius <- bw.stoyan(X)
  if(nsim > 1) {
    result <- vector(mode="list", length=nsim)
    for(isim in 1:nsim)
      result[[isim]] <- rjitter(X, radius=radius,
                                retry=retry, giveup=giveup, ...)
    names(result) <- paste("Simulation", 1:nsim)
    return(as.solist(result))
  }
  nX <- npoints(X)
  if(nX == 0) return(X)
  W <- X$window
  if(!retry) {
    ## points outside window are lost
    D <- runifdisc(nX, radius=radius)
    xnew <- X$x + D$x
    ynew <- X$y + D$y
    ok <- inside.owin(xnew, ynew, W)
    return(ppp(xnew[ok], ynew[ok], window=W))
  }
  ## retry = TRUE: condition on points being inside window
  undone <- rep.int(TRUE, nX)
  while(any(undone)) {
    giveup <- giveup - 1
    if(giveup <= 0)
      return(X)
    Y <- X[undone]
    D <- runifdisc(Y$n, radius=radius)
    xnew <- Y$x + D$x
    ynew <- Y$y + D$y
    ok <- inside.owin(xnew, ynew, W)
    if(any(ok)) {
      changed <- seq_len(nX)[undone][ok]
      X$x[changed] <- xnew[ok]
      X$y[changed] <- ynew[ok]
      undone[changed] <- FALSE
    }
  }
  return(if(drop) X else solist(X))
}

