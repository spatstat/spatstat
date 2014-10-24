#
#    random.R
#
#    Functions for generating random point patterns
#
#    $Revision: 4.63 $   $Date: 2014/10/24 00:22:30 $
#
#
#    runifpoint()      n i.i.d. uniform random points ("binomial process")
#
#    runifpoispp()     uniform Poisson point process
#
#    rpoispp()         general Poisson point process (thinning method)
#
#    rpoint()          n independent random points (rejection/pixel list)
#
#    rMaternI()        Mat'ern model I 
#    rMaternII()       Mat'ern model II
#    rSSI()            Simple Sequential Inhibition process
#
#    rNeymanScott()    Neyman-Scott process (generic)
#    rMatClust()       Mat'ern cluster process
#    rThomas()         Thomas process
#
#    rthin()           independent random thinning
#    rjitter()         random perturbation
#
#    Examples:
#          u01 <- owin(0:1,0:1)
#          plot(runifpoispp(100, u01))
#          X <- rpoispp(function(x,y) {100 * (1-x/2)}, 100, u01)
#          X <- rpoispp(function(x,y) {ifelse(x < 0.5, 100, 20)}, 100)
#          plot(X)
#          plot(rMaternI(100, 0.02))
#          plot(rMaternII(100, 0.05))
#

"runifrect" <-
  function(n, win=owin(c(0,1),c(0,1)))
{
  # no checking
      x <- runif(n, min=win$xrange[1], max=win$xrange[2])
      y <- runif(n, min=win$yrange[1], max=win$yrange[2])  
      return(ppp(x, y, window=win, check=FALSE))
}

"runifdisc" <-
  function(n, radius=1, centre=c(0,0), ...)
{
  # i.i.d. uniform points in the disc of radius r and centre (x,y)
  disque <- disc(centre=centre, radius=radius, ...)
  theta <- runif(n, min=0, max= 2 * pi)
  s <- sqrt(runif(n, min=0, max=radius^2))
  return(ppp(centre[1] + s * cos(theta),
             centre[2] + s * sin(theta),
             window=disque, check=FALSE))
}


"runifpoint" <-
  function(n, win=owin(c(0,1),c(0,1)), giveup=1000, warn=TRUE)
{
    win <- as.owin(win)
    
    check.1.integer(n)
    stopifnot(n >= 0)

    if(n == 0)
      return(ppp(numeric(0), numeric(0), window=win))

    if(warn) {
      nhuge <- spatstat.options("huge.npoints")
      if(n > nhuge)
        warning(paste("Attempting to generate", n, "random points"))
    }

    switch(win$type,
           rectangle = {
             return(runifrect(n, win))
           },
           mask = {
             dx <- win$xstep
             dy <- win$ystep
             # extract pixel coordinates and probabilities
             rxy <- rasterxy.mask(win, drop=TRUE)
             xpix <- rxy$x
             ypix <- rxy$y
             # select pixels with equal probability
             id <- sample(seq_along(xpix), n, replace=TRUE)
             # extract pixel centres and randomise within pixels
             x <- xpix[id] + runif(n, min= -dx/2, max=dx/2)
             y <- ypix[id] + runif(n, min= -dy/2, max=dy/2)
             return(ppp(x, y, window=win, check=FALSE))
           },
           polygonal={
             # rejection method
             # initialise empty pattern
             x <- numeric(0)
             y <- numeric(0)
             X <- ppp(x, y, window=win)
             #
             # rectangle in which trial points will be generated
             box <- boundingbox(win)
             # 
             ntries <- 0
             repeat {
               ntries <- ntries + 1
               # generate trial points in batches of n
               qq <- runifrect(n, box) 
               # retain those which are inside 'win'
               qq <- qq[win]
               # add them to result
               X <- superimpose(X, qq, W=win)
               # if we have enough points, exit
               if(X$n > n) 
                 return(X[1:n])
               else if(X$n == n)
                 return(X)
               # otherwise get bored eventually
               else if(ntries >= giveup)
                 stop(paste("Gave up after", giveup * n, "trials,",
                            X$n, "points accepted"))
             }
           })
    stop("Unrecognised window type")
}

"runifpoispp" <-
function(lambda, win = owin(c(0,1),c(0,1))) {
    win <- as.owin(win)
    if(!is.numeric(lambda) || length(lambda) > 1 ||
       !is.finite(lambda) || lambda < 0)
      stop("Intensity lambda must be a single finite number >= 0")

    if(lambda == 0) # return empty pattern
      return(ppp(numeric(0), numeric(0), window=win))

    # generate Poisson process in enclosing rectangle 
    box <- boundingbox(win)
    mean <- lambda * area(box)
    n <- rpois(1, mean)
    X <- runifpoint(n, box)

    # trim to window
    if(win$type != "rectangle")
      X <- X[win]  

    return(X)
}

rpoint <- function(n, f, fmax=NULL,
                   win=unit.square(), ..., giveup=1000,verbose=FALSE) {
  
  if(missing(f) || (is.numeric(f) && length(f) == 1))
    # uniform distribution
    return(runifpoint(n, win, giveup))
  
  # non-uniform distribution....
  
  if(!is.function(f) && !is.im(f))
    stop(paste(sQuote("f"),
               "must be either a function or an",
               sQuote("im"), "object"))
  
  if(is.im(f)) {
    # ------------ PIXEL IMAGE ---------------------
    wf <- as.owin(f)
    if(n == 0)
      return(ppp(numeric(0), numeric(0), window=wf))
    w <- as.mask(wf)
    M <- w$m
    dx <- w$xstep
    dy <- w$ystep
    # extract pixel coordinates and probabilities
    rxy <- rasterxy.mask(w, drop=TRUE)
    xpix <- rxy$x
    ypix <- rxy$y
    ppix <- as.vector(f$v[M]) # not normalised - OK
    # select pixels
    id <- sample(length(xpix), n, replace=TRUE, prob=ppix)
    # extract pixel centres and randomise within pixels
    x <- xpix[id] + runif(n, min= -dx/2, max=dx/2)
    y <- ypix[id] + runif(n, min= -dy/2, max=dy/2)
    return(ppp(x, y, window=wf, check=FALSE))
  }

  # ------------ FUNCTION  ---------------------  
  # Establish parameters for rejection method

  verifyclass(win, "owin")
  if(n == 0)
    return(ppp(numeric(0), numeric(0), window=win))
  
  if(is.null(fmax)) {
    # compute approx maximum value of f
    imag <- as.im(f, win, ...)
    summ <- summary(imag)
    fmax <- summ$max + 0.05 * diff(summ$range)
  }
  irregular <- (win$type != "rectangle")
  box <- boundingbox(win)
  X <- ppp(numeric(0), numeric(0), window=win)
  
  ntries <- 0

  # generate uniform random points in batches
  # and apply the rejection method.
  # Collect any points that are retained in X

  repeat{
    ntries <- ntries + 1
    # proposal points
    prop <- runifrect(n, box)
    if(irregular)
      prop <- prop[win]
    if(prop$n > 0) {
      fvalues <- f(prop$x, prop$y, ...)
      paccept <- fvalues/fmax
      u <- runif(prop$n)
      # accepted points
      Y <- prop[u < paccept]
      if(Y$n > 0) {
        # add to X
        X <- superimpose(X, Y, W=win)
        if(X$n >= n) {
          # we have enough!
          if(verbose)
            cat(paste("acceptance rate = ",
                      round(100 * X$n/(ntries * n), 2), "%\n"))
          return(X[1:n])
        }
      }
    }
    if(ntries > giveup)
      stop(paste("Gave up after",giveup * n,"trials with",
                 X$n, "points accepted"))
  }
  invisible(NULL)
}

"rpoispp" <-
  function(lambda, lmax=NULL, win = owin(c(0,1),c(0,1)), ...) {
    # arguments:
    #     lambda  intensity: constant, function(x,y,...) or image
    #     lmax     maximum possible value of lambda(x,y,...)
    #     win     default observation window (of class 'owin')
    #   ...       arguments passed to lambda(x, y, ...)

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
    
    if(is.numeric(lambda)) 
      # uniform Poisson
      return(runifpoispp(lambda, win))

    # inhomogeneous Poisson
    # perform thinning of uniform Poisson

    if(is.null(lmax)) {
      imag <- as.im(lambda, win, ...)
      summ <- summary(imag)
      lmax <- summ$max + 0.05 * diff(summ$range)
    } 

    if(is.function(lambda)) {
      X <- runifpoispp(lmax, win)  # includes sanity checks on `lmax'
      if(X$n == 0) return(X)
      prob <- lambda(X$x, X$y, ...)/lmax
      u <- runif(X$n)
      retain <- (u <= prob)
      X <- X[retain, ]
      return(X)
    }
    if(is.im(lambda)) {
      X <- runifpoispp(lmax, win)
      if(X$n == 0) return(X)
      prob <- lambda[X]/lmax
      u <- runif(X$n)
      retain <- (u <= prob)
      X <- X[retain, ]
      return(X)
    }
    stop(paste(sQuote("lambda"), "must be a constant, a function or an image"))
}
    
"rMaternI" <-
  function(kappa, r, win = owin(c(0,1),c(0,1)), stationary=TRUE)
{
  win <- as.owin(win)
  stopifnot(is.numeric(r) && length(r) == 1)
  if(stationary) {
    # generate in a larger window
    bigbox <- grow.rectangle(as.rectangle(win), r)
    X <- rpoispp(kappa, win=bigbox)
  } else {
    X <- rpoispp(kappa, win=win)
  }
  if(npoints(X) > 1) {
    d <- nndist(X)
    X <- X[d > r]
  }
  if(stationary)
    X <- X[win]
  return(X)
}
    
"rMaternII" <-
  function(kappa, r, win = owin(c(0,1),c(0,1)), stationary=TRUE)
{
  win <- as.owin(win)
  stopifnot(is.numeric(r) && length(r) == 1)
  if(stationary) {
    bigbox <- grow.rectangle(as.rectangle(win), r)
    X <- rpoispp(kappa, win=bigbox)
  } else {
    X <- rpoispp(kappa, win=win)
  }

  nX <- npoints(X)
  if(nX > 1) {
    # matrix of squared pairwise distances
    d2 <- pairdist(X, squared=TRUE)
    close <- (d2 <= r^2)
    # random order 1:n
    age <- sample(seq_len(nX), nX, replace=FALSE)
    earlier <- outer(age, age, ">")
    conflict <- close & earlier
    # delete <- apply(conflict, 1, any)
    delete <- matrowany(conflict)
    X <- X[!delete]
  }
  if(stationary)
    X <- X[win]
  return(X)
}
  
"rSSI" <-
  function(r, n=Inf, win = square(1), 
           giveup = 1000, x.init=NULL)
{
  win.given <- !missing(win) && !is.null(win)
  stopifnot(is.numeric(r) && length(r) == 1 && r >= 0)
  stopifnot(is.numeric(n) && length(n) == 1 && n >= 0)
     # Simple Sequential Inhibition process
     # fixed number of points
     # Naive implementation, proposals are uniform
  if(is.null(x.init)) {
    # start with empty pattern in specified window
    win <- as.owin(win)
    x.init <- ppp(numeric(0),numeric(0), window=win)
  } else {
    # start with specified pattern
    stopifnot(is.ppp(x.init))
    if(!win.given) {
      win <- as.owin(x.init)
    } else {
      # check compatibility of windows
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
      return(x.init)
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
    qq <- runifpoint(1, win)
    x <- qq$x[1]
    y <- qq$y[1]
    if(X$n == 0 || all(((x - X$x)^2 + (y - X$y)^2) > r2))
      X <- superimpose(X, qq, W=win)
    if(X$n >= n)
      return(X)
  }
  if(!is.infinite(n))
    warning(paste("Gave up after", giveup,
                  "attempts with only", X$n, "points placed out of", n))
  return(X)
}

"rPoissonCluster" <-
  function(kappa, rmax, rcluster, win = owin(c(0,1),c(0,1)), ..., lmax=NULL)
{
  # Generic Poisson cluster process
  # Implementation for bounded cluster radius
  #
  # 'rcluster' is a function(x,y) that takes the coordinates
  # (x,y) of the parent point and generates a list(x,y) of offspring
  #
  # "..." are arguments to be passed to 'rcluster()'
  #

  win <- as.owin(win)
  
  # Generate parents in dilated window
  frame <- boundingbox(win)
  dilated <- owin(frame$xrange + c(-rmax, rmax),
                  frame$yrange + c(-rmax, rmax))
  if(is.im(kappa) && !is.subset.owin(dilated, as.owin(kappa)))
    stop(paste("The window in which the image",
               sQuote("kappa"),
               "is defined\n",
               "is not large enough to contain the dilation of the window",
               sQuote("win")))
  parents <- rpoispp(kappa, lmax=lmax, win=dilated)
  #
  result <- NULL
  # generate clusters
  np <- parents$n
  if(np > 0) {
    xparent <- parents$x
    yparent <- parents$y
    for(i in seq_len(np)) {
      # generate random offspring of i-th parent point
      cluster <- rcluster(xparent[i], yparent[i], ...)
      if(!inherits(cluster, "ppp"))
        cluster <- ppp(cluster$x, cluster$y, window=frame, check=FALSE)
      # skip if cluster is empty
      if(cluster$n > 0) {
        # trim to window
        cluster <- cluster[win]
        if(is.null(result)) {
          # initialise offspring pattern and offspring-to-parent map
          result <- cluster
          parentid <- rep.int(1, cluster$n)
        } else {
          # add to pattern
          result <- superimpose(result, cluster, W=win)
          # update offspring-to-parent map
          parentid <- c(parentid, rep.int(i, cluster$n))
        }
      }
    }
  } else {
    # no parents - empty pattern
    result <- ppp(numeric(0), numeric(0), window=win)
    parentid <- integer(0)
  }

  attr(result, "parents") <- parents
  attr(result, "parentid") <- parentid
  
  return(result)
}  

rGaussPoisson <-
  function(kappa, r, p2, win=owin(c(0,1), c(0,1)))
{
  # Gauss-Poisson process
  oneortwo <- function(x0, y0, radius, p2) {
    if(runif(1) > p2) 
      # one point
      return(list(x=x0, y=y0))
    # two points
    theta <- runif(1, min=0, max=2*pi)
    return(list(x=x0+c(-1,1)*radius*cos(theta),
                y=y0+c(-1,1)*radius*sin(theta)))
  }
  result <- rPoissonCluster(kappa, 1.05 * r, oneortwo,
                            win, radius=r/2, p2=p2)
  return(result)
  
}

  
rstrat <- function(win=square(1), nx, ny=nx, k=1) {
  win <- as.owin(win)
  stopifnot(nx >= 1 && ny >= 1)
  stopifnot(k >= 1)
  xy <- stratrand(win, nx, ny, k)
  Xbox <- ppp(xy$x, xy$y, win$xrange, win$yrange, check=FALSE)
  X <- Xbox[win]
  return(X)
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
    nx <- length(x0)
  } else stop("Need either nx or dx")
  # determine y grid
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
    ny <- length(y0)
  }
  return(list(x0=x0, y0=y0, nx=nx, ny=ny, dx=dx, dy=dy))
}
  
rsyst <- function(win=square(1), nx=NULL, ny=nx, ..., dx=NULL, dy=dx) {
  win <- as.owin(win)
  xr <- win$xrange
  yr <- win$yrange
  # determine grid coordinates 
  if(missing(ny)) ny <- NULL
  if(missing(dy)) dy <- NULL
  g <- xy.grid(xr, yr, nx, ny, dx, dy)
  x0 <- g$x0
  y0 <- g$y0
  dx <- g$dx
  dy <- g$dy
  # assemble grid and randomise location
  xy0 <- expand.grid(x=x0, y=y0)
  x <- xy0$x + runif(1, min = 0, max = dx)
  y <- xy0$y + runif(1, min = 0, max = dy)
  Xbox <- ppp(x, y, xr, yr, check=FALSE)
  # trim to window
  X <- Xbox[win]
  return(X)
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
  k <- ifelse(u < 1/N, 0, ifelse(u < (1 - pN), 1, N))
  return(k)
}

rcell <- function(win=square(1), nx=NULL, ny=nx, ..., dx=NULL, dy=dx, N=10) {
  win <- as.owin(win)
  xr <- win$xrange
  yr <- win$yrange
  # determine grid coordinates 
  if(missing(ny)) ny <- NULL
  if(missing(dy)) dy <- NULL
  g <- xy.grid(xr, yr, nx, ny, dx, dy)
  nx <- g$nx
  ny <- g$ny
  x0 <- g$x0
  y0 <- g$y0
  dx <- g$dx
  dy <- g$dy
  # generate pattern
  x <- numeric(0)
  y <- numeric(0)
  for(ix in seq_len(nx))
    for(iy in seq_len(ny)) {
      nij <- rcellnumber(1, N)
      x <- c(x, x0[ix] + runif(nij, min=0, max=dx))
      y <- c(y, y0[iy] + runif(nij, min=0, max=dy))
    }
  Xbox <- ppp(x, y, xr, yr, check=FALSE)
  X <- Xbox[win]
  return(X)
}


rthin <- function(X, P, ...) {
  verifyclass(X, "ppp")

  nX <- npoints(X)
  if(nX == 0) return(X)

  if(is.numeric(P)) {
    # vector of retention probabilities
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
    # function - evaluate it at points of X
    pX <- P(X$x, X$y, ...)
    if(length(pX) != nX)
      stop("Function P returned a vector of incorrect length")
    if(!is.numeric(pX))
      stop("Function P returned non-numeric values")
    if(any(is.na(pX)))
      stop("Function P returned some NA values")
    prange <- range(pX)
  } else if(is.im(P)) {
    # image - look it up
    if(!(P$type %in% c("integer", "real")))
      stop("Values of image P should be numeric")
    pX <- P[X, drop=FALSE]
    if(any(is.na(pX)))
      stop("some points of X lie outside the domain of image P")
  } else
  stop("Unrecognised format for P")

  if(min(pX) < 0) stop("some probabilities are negative")
  if(max(pX) > 1) stop("some probabilities are greater than 1")

  retain <- (runif(length(pX)) < pX)

  Y <- X[retain]
  
  # also handle offspring-to-parent map if present
  if(!is.null(parentid <- attr(X, "parentid")))
    attr(Y, "parentid") <- parentid[retain]
  
  return(Y)
}


# rjitter

rjitter <- function(X, radius, retry=TRUE, giveup=10000) {
  verifyclass(X, "ppp")
  if(missing(radius) || is.null(radius))
    radius <- bw.stoyan(X)
  nX <- npoints(X)
  if(nX == 0) return(X)
  W <- X$window
  if(!retry) {
    # points outside window are lost
    D <- runifdisc(nX, radius=radius)
    xnew <- X$x + D$x
    ynew <- X$y + D$y
    ok <- inside.owin(xnew, ynew, W)
    return(ppp(xnew[ok], ynew[ok], window=W))
  }
  # retry = TRUE: condition on points being inside window
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
  return(X)
}

