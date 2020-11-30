#'
#'   randombasic.R
#'
#'   Basic random generators, needed in spatstat.geom
#' 
#'    runifrect()       special case of rectangle
#'    rsyst()           systematic random (randomly-displaced grid)
#'    rjitter()         random perturbation
#'
#'   $Revision: 1.5 $  $Date: 2020/11/30 08:05:10 $


simulationresult <- function(resultlist, nsim, drop, NameBase="Simulation") {
  if(nsim == 1 && drop)
    return(resultlist[[1L]])
  #' return 'solist' if appropriate, otherwise 'anylist'
  return(as.solist(resultlist, .NameBase=NameBase, demote=TRUE))
}

runifrect <- function(n, win=owin(c(0,1),c(0,1)), nsim=1, drop=TRUE)
{
  ## no checking
  xr <- win$xrange
  yr <- win$yrange
  result <- vector(mode="list", length=nsim)
  for(isim in 1:nsim) {
    x <- runif(n, min=xr[1], max=xr[2])
    y <- runif(n, min=yr[1], max=yr[2])
    result[[isim]] <- ppp(x, y, window=win, check=FALSE)
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}

rsyst <- function(win=square(1), nx=NULL, ny=nx, ..., dx=NULL, dy=dx,
                  nsim=1, drop=TRUE) {
  if(!missing(nsim)) {
    check.1.integer(nsim)
    stopifnot(nsim >= 1)
  }
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
  result <- vector(mode="list", length=nsim)
  for(isim in 1:nsim) {
    x <- xy0$x + runif(1, min = 0, max = dx)
    y <- xy0$y + runif(1, min = 0, max = dy)
    Xbox <- ppp(x, y, xr, yr, check=FALSE)
    ## trim to window
    result[[isim]] <- Xbox[win]
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
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

## rjitter

rjitter <- function(X, radius, retry=TRUE, giveup=10000, ...,
                    nsim=1, drop=TRUE) {
  verifyclass(X, "ppp")
  check.1.integer(nsim)
  stopifnot(nsim >= 1)
  nX <- npoints(X)
  if(nX == 0) {
    result <- rep(list(X), nsim)
    result <- simulationresult(result, nsim, drop)
    return(result)
  }
  W <- X$window
  if(missing(radius) || is.null(radius)) {
    ## Stoyan rule of thumb
    bws <- 0.15/sqrt(5 * nX/area(W))
    radius <- min(bws, shortside(Frame(W)))
  }
  
  result <- vector(mode="list", length=nsim)
  for(isim in 1:nsim) {
    if(!retry) {
      ## points outside window are lost
      rD <- sqrt(runif(nX, max=radius^2))
      aD <- runif(nX, max= 2 * pi)
      xnew <- X$x + rD * cos(aD)
      ynew <- X$y + rD * sin(aD)
      ok <- inside.owin(xnew, ynew, W)
      result[[isim]] <- ppp(xnew[ok], ynew[ok], window=W, check=FALSE)
    } else {
      ## retry = TRUE: condition on points being inside window
      undone <- rep.int(TRUE, nX)
      triesleft <- giveup
      Xshift <- X
      while(any(undone)) {
        triesleft <- triesleft - 1
        if(triesleft <= 0) 
	  break
        Y <- Xshift[undone]
        nY <- npoints(Y)
        rD <- sqrt(runif(nY, max=radius^2))
        aD <- runif(nY, max= 2 * pi)
        xnew <- Y$x + rD * cos(aD)
        ynew <- Y$y + rD * sin(aD)
        ok <- inside.owin(xnew, ynew, W)
        if(any(ok)) {
          changed <- which(undone)[ok]
          Xshift$x[changed] <- xnew[ok]
          Xshift$y[changed] <- ynew[ok]
          undone[changed] <- FALSE
        }
      }
      result[[isim]] <- Xshift
    }
  }
  result <- simulationresult(result, nsim, drop)
  return(result)
}

