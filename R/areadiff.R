#
# areadiff.R
#
#  $Revision: 1.28 $  $Date: 2013/11/02 01:53:09 $
#
# Computes sufficient statistic for area-interaction process
#
# Invokes areadiff.c
#
# areaLoss = area lost by removing X[i] from X

areaLoss <- function(X, r, ..., W=as.owin(X),
                     subset=NULL, exact=FALSE,
                     ngrid=spatstat.options("ngrid.disc")) {
  if(exact)
    areaLoss.diri(X, r, ..., W=W, subset=subset)
  else
    areaLoss.grid(X, r, ..., W=W, subset=subset, ngrid=ngrid)
}

# areaGain = area gained by adding u[i] to X

areaGain <- function(u, X, r, ..., W=as.owin(X), exact=FALSE,
                     ngrid=spatstat.options("ngrid.disc")) {
  if(exact)
    areaGain.diri(u, X, r, ..., W=W)
  else
    areaGain.grid(u, X, r, W=W, ngrid=ngrid)
}


#////////////////////////////////////////////////////////////
#    algorithms using Dirichlet tessellation
#///////////////////////////////////////////////////////////

areaLoss.diri <- function(X, r, ..., W=as.owin(X), subset=NULL) {
  stopifnot(is.ppp(X))
  npts <- npoints(X)
  if(is.matrix(r)) {
    if(sum(dim(r) > 1) > 1)
      stop("r should be a vector or single value")
    r <- as.vector(r)
  }
  nr <- length(r)
  if(npts == 0)
    return(matrix(, nrow=0, ncol=nr))
  else if(npts == 1) 
    return(matrix(discpartarea(X, r, W), nrow=1))
  # set up output array
  indices <- 1:npts
  if(!is.null(subset))
    indices <- indices[subset]
  out <- matrix(, nrow=length(indices), ncol=nr)
  #
  w <- X$window
  pir2 <- pi * r^2
  # dirichlet neighbour relation in entire pattern 
  dd <- deldir(X$x, X$y, rw=c(w$xrange, w$yrange))
  a <- dd$delsgs[,5]
  b <- dd$delsgs[,6]
  for(k in seq_along(indices)) {
    i <- indices[k]
    # find all Delaunay neighbours of i 
    jj <- c(b[a==i], a[b==i])
    jj <- sort(unique(jj))
    # extract only these points
    Yminus <- X[jj]
    Yplus  <- X[c(jj, i)]
    # dilate
    aplus <- dilated.areas(Yplus, r, W, exact=TRUE)
    aminus <- dilated.areas(Yminus, r, W, exact=TRUE)
    areas <- aplus - aminus
    # area/(pi * r^2) must be positive and nonincreasing
    y <- ifelseAX(r == 0, 1, areas/pir2)
    y <- pmin.int(1, y)
    ok <- is.finite(y)
    y[ok] <- rev(cummax(rev(y[ok])))
    areas <- pmax.int(0, y * pir2)
    # save
    out[k, ] <- areas
  }
  return(out)
}

areaGain.diri <- function(u, X, r, ..., W=as.owin(X)) {
  stopifnot(is.ppp(X))
  Y <- as.ppp(u, W=W)
  nX <- X$n
  nY <- Y$n
  if(is.matrix(r)) {
    if(sum(dim(r) > 1) > 1)
      stop("r should be a vector or single value")
    r <- as.vector(r)
  }
  nr <- length(r)
  if(nY == 0)
    return(matrix(, nrow=0, ncol=nr))
  if(nX == 0)
    return(matrix(pi * r^2, nrow=nY, ncol=nr, byrow=TRUE))
  cat(paste("areaGain,", nY, "points,", nr, "r values\n"))
  out <- matrix(0, nrow=nY, ncol=nr)
  pir2 <- pi * r^2
  wbox <- as.rectangle(as.owin(X))
  #
  for(i in 1:nY) {
    progressreport(i, nY)
    V <- superimpose(Y[i], X, W=wbox, check=FALSE)
    # Dirichlet neighbour relation for V
    dd <- deldir(V$x, V$y, rw=c(wbox$xrange, wbox$yrange))
    aa <- dd$delsgs[,5]
    bb <- dd$delsgs[,6]
    # find all Delaunay neighbours of Y[1] in V
    jj <- c(bb[aa==1], aa[bb==1])
    jj <- sort(unique(jj))
    # extract only these points
    Zminus <- V[jj]
    Zplus  <- V[c(1, jj)]
    # dilate
    aplus <- dilated.areas(Zplus, r, W, exact=TRUE)
    aminus <- dilated.areas(Zminus, r, W, exact=TRUE)
    areas <- aplus - aminus
    # area/(pi * r^2) must be in [0,1] and nonincreasing
    y <- ifelseAX(r == 0, 1, areas/pir2)
    y <- pmin.int(1, y)
    ok <- is.finite(y)
    y[ok] <- rev(cummax(rev(y[ok])))
    areas <- pmax.int(0, y * pir2)
    # save
    out[i,] <- areas
  }
  return(out)
}

#////////////////////////////////////////////////////////////////////////
#    alternative implementations using grid counting in C
#////////////////////////////////////////////////////////////////////////

areaGain.grid <- function(u, X, r, ..., W=NULL, ngrid=spatstat.options("ngrid.disc")) {
  verifyclass(X, "ppp")
  u <- as.ppp(u, W=as.owin(X))
  stopifnot(is.numeric(r) && all(is.finite(r)) && all(r >= 0))
  #
  nu <- u$n
  nr <- length(r)
  if(nr == 0)
    return(numeric(0))
  rmax <- max(r)
  #
  constrain <- !is.null(W)
  if(constrain && (W$type != "rectangle")) {
    # Constrained to an irregular window
    # initialise to value for small-r
    result <- matrix(pi * r^2, nrow=nu, ncol=nr, byrow=TRUE)    
    # vector of radii below which b(u,r) is disjoint from U(X,r)
    rcrit.u <- nncross(u, X, what="dist")/2
    rcrit.min <- min(rcrit.u)
    # Use distance transform and set covariance
    D <- distmap(X, ...)
    DW <- D[W, drop=FALSE]
    # distance from (0,0) - thresholded to make digital discs
    discWin <- owin(c(-rmax,rmax),c(-rmax,rmax))
    discWin <- as.mask(discWin, eps=min(D$xstep, rmax/4))
    rad <- as.im(function(x,y){sqrt(x^2+y^2)}, W=discWin)
    # 
    for(j in which(r > rcrit.min)) {
      # rj is above the critical radius rcrit.u[i] for at least one point u[i]
      rj <- r[j]
      if(any(above <- (rj > rcrit.u))) {
        Uncovered  <- levelset(DW, rj, ">")
        DiscRj     <- levelset(rad, rj, "<=")
        AreaGainIm <- setcov(Uncovered, DiscRj)
        result[above, j] <- safelookup(AreaGainIm, u[above])
      }
    }
    return(result)
  }
  #
  #
  xx <- X$x
  yy <- X$y
  result <- matrix(, nrow=nu, ncol=nr)
  DUP <- spatstat.options("dupC")
  #
  for(i in 1:nu) {
    # shift u[i] to origin
    xu <- u$x[i]
    yu <- u$y[i]
    xshift <- xx - xu
    yshift <- yy - yu
    # find points within distance 2 rmax of origin
    close <- (xshift^2 + yshift^2 < 4 * rmax^2)
    nclose <- sum(close)
    # invoke C routine
    if(!constrain) {
      z <- .C("areadifs",
              rad = as.double(r),
              nrads = as.integer(nr),
              x   = as.double(xshift[close]),
              y   = as.double(yshift[close]),
              nn  = as.integer(nclose),
              ngrid = as.integer(ngrid),
              answer = as.double(numeric(nr)),
              DUP=DUP)
#              PACKAGE="spatstat")
      result[i,] <- z$answer
    } else {
      z <- .C("areaBdif",
              rad = as.double(r),
              nrads = as.integer(nr),
              x   = as.double(xshift[close]),
              y   = as.double(yshift[close]),
              nn  = as.integer(nclose),
              ngrid = as.integer(ngrid),
              x0 = as.double(W$xrange[1] - xu),
              y0 = as.double(W$yrange[1] - yu),
              x1 = as.double(W$xrange[2] - xu),
              y1 = as.double(W$yrange[2] - yu),
              answer = as.double(numeric(nr)),
              DUP=DUP)
#              PACKAGE="spatstat")
      result[i,] <- z$answer
    }
  }
  return(result)
}

areaLoss.grid <- function(X, r, ...,
                          W=as.owin(X), subset=NULL,
                          method = c("count", "distmap"),
                          ngrid = spatstat.options("ngrid.disc"),
                          exact = FALSE) {
  verifyclass(X, "ppp")
  n <- npoints(X)
  nr <- length(r)
  indices <- if(is.null(subset)) 1:n else (1:n)[subset]
  answer <- matrix(, nrow=length(indices), ncol=nr)
  if(missing(method)) {
    method <- if(nr <= 20 || exact) "count" else "distmap"
  } else method <- match.arg(method)
  switch(method,
         count = {
           # one value of r: use grid-counting
           for(k in seq_along(indices)) {
             i <- indices[k]
             answer[k,] <- areaGain(X[i], X[-i], r, W=W,
                                    ngrid=ngrid, exact=exact)
           }
         },
         distmap = {
           # Many values of r: use distance transform
           D <- distmap(X, ...)
           DW <- D[W, drop=FALSE]
           a <- area(Window(DW))
           # empirical cdf of distance values
           FW <- ecdf(DW[drop=TRUE])
           # radii below which there are no overlaps
           rcrit <- nndist(X)/2
           for(k in seq_along(indices)) {
             i <- indices[k]
             Di <- distmap(X[-i], ...)
             FiW <- ecdf(Di[W, drop=TRUE])
             answer[k, ] <-
               ifelseXY(r > rcrit[i], a * (FW(r) - FiW(r)), pi * r^2)
           }
         })
  return(answer)
}

