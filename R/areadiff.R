#
# areadiff.R
#
#  $Revision: 1.38 $  $Date: 2019/01/20 08:46:55 $
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
    areaGain.grid(u, X, r, W=W, ..., ngrid=ngrid)
}

#////////////////////////////////////////////////////////////
#    algorithms using polygon geometry
#///////////////////////////////////////////////////////////

areaLoss.poly <- function(X, r, ..., W=as.owin(X), subset=NULL,
                            splitem=TRUE) {
  check.1.real(r)
  nX <- npoints(X)
  if(r <= 0 || nX == 0) return(numeric(nX))
  cooX <- coords(X)
  if(useW <- is.owin(W))
    W <- as.polygonal(W)
  #' initialise result
  result <- rep(pi * r^2, nX)
  wanted <- 1:nX
  if(!is.null(subset))
    wanted <- wanted[subset]
  #' split into connected components
  if(splitem) {
    Y <- connected(X, 2 * r)
    Z <- split(Y)
    V <- lapply(Z, areaLoss.poly, r=r, W=W, splitem=FALSE)
    return(unsplit(V, marks(Y))[wanted])
  }
  #' determine which pairs of points interact
  cl <- closepairs(X, 2 * r, what="indices")
  if(length(cl$i) == 0)
    return(result[wanted])
  #' determine scale parameters for polyclip
  p <- commonPolyclipArgs(Frame(X))
  #' template disc
  ball0 <- disc(r, c(0,0), ...)
  #' discs centred on data points
  balls <- vector(mode="list", length=nX)
  for(i in seq_len(nX))
    balls[[i]] <- shift(ball0, vec=cooX[i,])
  balls <- as.solist(balls, check=FALSE)
  #' start computin'
  for(i in wanted) {
    jj <- cl$j[cl$i == i]
    nn <- length(jj)
    if(nn > 0) {
      #' union of balls close to i
      u <- if(nn == 1) balls[[ jj ]] else union.owin(balls[jj], p=p)
      #' subtract from ball i
      v <- setminus.owin(balls[[i]], u)
      #' clip to window
      if(useW) 
        v <- intersect.owin(v, W)
      #' compute
      result[i] <- area(v)
    }
  }
  return(result[wanted])
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
  #' set up output array
  indices <- 1L:npts
  if(!is.null(subset))
    indices <- indices[subset]
  out <- matrix(0, nrow=length(indices), ncol=nr)
  #' handle duplicate points
  retain <- !duplicated(X)
  getzero <- (multiplicity(X) > 1)
  uX <- X[retain]
  newserial <- cumsum(retain)
  # dirichlet neighbour relation in entire pattern 
  w <- X$window
  dd <- deldir(uX$x, uX$y, rw=c(w$xrange, w$yrange))
  a <- dd$delsgs[,5L]
  b <- dd$delsgs[,6L]
  pir2 <- pi * r^2
  for(k in seq_along(indices)) {
    ind <- indices[k]
    if(!getzero[ind]) {
      #' find serial number in uX
      i <- newserial[ind]
      #' find all Delaunay neighbours of i 
      jj <- c(b[a==i], a[b==i])
      jj <- sortunique(jj)
      #' extract only these points
      Yminus <- uX[jj]
      Yplus  <- uX[c(jj, i)]
      #' dilate
      aplus <- dilated.areas(Yplus, r, W, exact=TRUE, ...)
      aminus <- dilated.areas(Yminus, r, W, exact=TRUE, ...)
      areas <- aplus - aminus
      #' area/(pi * r^2) must be positive and nonincreasing
      y <- ifelseAX(r == 0, 1, areas/pir2)
      y <- pmin.int(1, y)
      ok <- is.finite(y)
      y[ok] <- rev(cummax(rev(y[ok])))
      areas <- pmax.int(0, y * pir2)
      #' save
      out[k, ] <- areas
    }
  }
  return(out)
}

areaGain.diri <- function(u, X, r, ..., W=as.owin(X), verbose=FALSE) {
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
  if(verbose)
    splat("areaGain,",
          nY, ngettext(nY, "point,", "points,"),
          nr, ngettext(nr, "rvalue", "r values"))
  out <- matrix(0, nrow=nY, ncol=nr)
  pir2 <- pi * r^2
  wbox <- as.rectangle(as.owin(X))
  #
  state <- list()
  for(i in 1L:nY) {
    if(verbose) state <- progressreport(i, nY, state=state)
    V <- superimpose(Y[i], X, W=wbox, check=FALSE)
    # Dirichlet neighbour relation for V
    dd <- deldir(V$x, V$y, rw=c(wbox$xrange, wbox$yrange))
    aa <- dd$delsgs[,5L]
    bb <- dd$delsgs[,6L]
    # find all Delaunay neighbours of Y[1] in V
    jj <- c(bb[aa==1L], aa[bb==1L])
    jj <- sortunique(jj)
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
    #' determine pixel resolution
    eps <- unclass(as.mask(Window(X), ...))[c("xstep", "ystep")]
    eps <- as.numeric(eps)
    eps <- eps * min(1, (rmax/4)/max(eps))
    #' Use distance transform and set covariance
    D <- distmap(X, eps=eps)
    DW <- D[W, drop=FALSE]
    # distance from (0,0) - thresholded to make digital discs
    discWin <- owin(c(-rmax,rmax),c(-rmax,rmax))
    discWin <- as.mask(discWin, eps=eps)
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
  #
  for(i in 1L:nu) {
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
              PACKAGE = "spatstat")
      result[i,] <- z$answer
    } else {
      z <- .C("areaBdif",
              rad = as.double(r),
              nrads = as.integer(nr),
              x   = as.double(xshift[close]),
              y   = as.double(yshift[close]),
              nn  = as.integer(nclose),
              ngrid = as.integer(ngrid),
              x0 = as.double(W$xrange[1L] - xu),
              y0 = as.double(W$yrange[1L] - yu),
              x1 = as.double(W$xrange[2L] - xu),
              y1 = as.double(W$yrange[2L] - yu),
              answer = as.double(numeric(nr)),
              PACKAGE = "spatstat")
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
  indices <- if(is.null(subset)) 1L:n else (1L:n)[subset]
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
                                    ngrid=ngrid, exact=exact, ...)
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

