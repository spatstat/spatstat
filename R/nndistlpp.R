#
# nndistlpp.R
#
#  $Revision: 1.3 $ $Date: 2013/10/21 02:35:05 $
#
# Methods for nndist, nnwhich, nncross for linear networks
#
# nndist.lpp
#   Calculates the nearest neighbour distances in the shortest-path metric
#   for a point pattern on a linear network.

nndist.lpp <- function(X, ..., method="C") {
  stopifnot(inherits(X, "lpp"))
  stopifnot(method %in% c("C", "interpreted"))
  #
  L <- X$domain
  Y <- as.ppp(X)
  n <- npoints(Y)
  #
  Lseg  <- L$lines
  Lvert <- L$vertices
  from  <- L$from
  to    <- L$to
  dpath <- L$dpath

  if(n == 0) return(numeric(0))
  if(n == 1) return(NA)
  
  # find nearest segment for each point
  # This is given by local coordinates, if available (spatstat >= 1.28-0)
  loco <- coords(X, local=TRUE, spatial=FALSE, temporal=FALSE)
  pro <- if(!is.null(seg <- loco$seg)) seg else nearestsegment(X, Lseg)

  if(method == "interpreted") {
    D <- pairdist(X, method="interpreted")
    diag(D) <- Inf
    return(apply(D, 1, min))
  } else {
    # C code
    # convert indices to start at 0
    from0 <- from - 1L
    to0   <- to - 1L
    segmap <- pro - 1L
    # upper bound on interpoint distance
    huge <- max(dpath) + 2 * max(lengths.psp(Lseg))
    # space for result
    ans <- double(n)
    # call C
    zz <- .C("linnndist",
             np = as.integer(n),
             xp = as.double(Y$x),
             yp = as.double(Y$y),
             nv = as.integer(Lvert$n),
             xv = as.double(Lvert$x),
             yv = as.double(Lvert$y),
             ns = as.double(L$n),
             from = as.integer(from0),
             to = as.integer(to0),
             dpath = as.double(dpath),
             segmap = as.integer(segmap),
             huge = as.double(huge),
             answer = as.double(ans))
    ans <- zz$answer
  }
  return(ans)
}

# nnwhich.lpp
# Identifies the nearest neighbours in the shortest-path metric
# for a point pattern on a linear network.
#

nnwhich.lpp <- function(X, ..., method="C") {
  stopifnot(inherits(X, "lpp"))
  stopifnot(method %in% c("C", "interpreted"))
  #
  L <- X$domain
  Y <- as.ppp(X)
  n <- npoints(Y)
  #
  Lseg  <- L$lines
  Lvert <- L$vertices
  from  <- L$from
  to    <- L$to
  dpath <- L$dpath
  
  if(n == 0) return(integer(0))
  if(n == 1) return(as.integer(NA))
  
  # find nearest segment for each point
  # This is given by local coordinates, if available (spatstat >= 1.28-0)
  loco <- coords(X, local=TRUE, spatial=FALSE, temporal=FALSE)
  pro <- if(!is.null(seg <- loco$seg)) seg else nearestsegment(X, Lseg)

  if(method == "interpreted") {
    D <- pairdist(X, method="interpreted")
    diag(D) <- Inf
    return(apply(D, 1, which.min))
  } else {
    # C code
    # convert indices to start at 0
    from0 <- from - 1
    to0   <- to - 1
    segmap <- pro - 1
    # upper bound on interpoint distance
    huge <- max(dpath) + 2 * max(lengths.psp(Lseg))
    # space for result
    nnd <- double(n)
    nnw <- integer(n)
    # call C
    zz <- .C("linnnwhich",
             np = as.integer(n),
             xp = as.double(Y$x),
             yp = as.double(Y$y),
             nv = as.integer(Lvert$n),
             xv = as.double(Lvert$x),
             yv = as.double(Lvert$y),
             ns = as.double(L$n),
             from = as.integer(from0),
             to = as.integer(to0),
             dpath = as.double(dpath),
             segmap = as.integer(segmap),
             huge = as.double(huge),
             nndist = as.double(nnd),
             nnwhich = as.integer(nnw))
    # convert C indexing to R indexing
    nnw <- zz$nnwhich + 1L
    # any zeroes occur if points have no neighbours.
    nnw[nnw == 0] <- NA
  }
  return(nnw)
}

# nncross.lpp
# Identifies the nearest neighbours in the shortest-path metric
# from one point pattern on a linear network to ANOTHER pattern
# on the SAME network.
#

nncross.lpp <- function(X, Y, iX=NULL, iY=NULL, what = c("dist", "which"), ..., method="C") {
  stopifnot(inherits(X, "lpp"))
  stopifnot(inherits(Y, "lpp"))
  what   <- match.arg(what, choices=c("dist", "which"), several.ok=TRUE)
  stopifnot(method %in% c("C", "interpreted"))
  check <- resolve.defaults(list(...), list(check=TRUE))$check
  #
  L <- as.linnet(X)
  if(check && !identical(L, as.linnet(Y))) 
    stop("X and Y are on different linear networks")
  #
  nX <- npoints(X)
  nY <- npoints(Y)
  P <- as.ppp(X)
  Q <- as.ppp(Y)
  #
  Lseg  <- L$lines
  Lvert <- L$vertices
  from  <- L$from
  to    <- L$to
  dpath <- L$dpath
  
  # deal with null cases
  if(nX == 0)
    return(data.frame(dist=numeric(0), which=integer(0))[, what])
  if(nY == 0)
    return(data.frame(dist=rep(Inf, nX), which=rep(NA_integer_, nX))[, what])

  # find nearest segment for each point
  Xpro <- coords(X, local=TRUE, spatial=FALSE, temporal=FALSE)$seg
  Ypro <- coords(Y, local=TRUE, spatial=FALSE, temporal=FALSE)$seg

  # handle serial numbers
  if(is.null(iX) != is.null(iY))
    stop("If one of iX, iY is given, then both must be given")
  exclude <- (!is.null(iX) || !is.null(iY))
  if(exclude) {
    stopifnot(is.integer(iX) && is.integer(iY))
    if(length(iX) != nX)
      stop("length of iX does not match the number of points in X")
    if(length(iY) != nY)
      stop("length of iY does not match the number of points in Y")
  }

  if(method == "interpreted") {
    D <- crossdist(X, method="interpreted")
    if(exclude)
      D[outer(iX, iY, "==")] <- Inf
    nnd <- if("dist" %in% what) apply(D, 1, min) else NA
    nnw <- if("which" %in% what) apply(D, 1, which.min) else NA
  } else {
    # C code
    # convert indices to start at 0
    from0 <- from - 1
    to0   <- to - 1
    Xsegmap <- Xpro - 1
    Ysegmap <- Ypro - 1
    # upper bound on interpoint distance
    huge <- max(dpath) + 2 * diameter(as.rectangle(as.owin(L)))
    # space for result
    nnd <- double(nX)
    nnw <- integer(nX)
    # call C
    if(!exclude) {
      zz <- .C("linndcross",
               np = as.integer(nX),
               xp = as.double(P$x),
               yp = as.double(P$y),
               nq = as.integer(nY),
               xq = as.double(Q$x),
               yq = as.double(Q$y),
               nv = as.integer(Lvert$n),
               xv = as.double(Lvert$x),
               yv = as.double(Lvert$y),
               ns = as.double(L$n),
               from = as.integer(from0),
               to = as.integer(to0),
               dpath = as.double(dpath),
               psegmap = as.integer(Xsegmap),
               qsegmap = as.integer(Ysegmap),
               huge = as.double(huge),
               nndist = as.double(nnd),
               nnwhich = as.integer(nnw))
    } else {
      zz <- .C("linndxcross",
               np = as.integer(nX),
               xp = as.double(P$x),
               yp = as.double(P$y),
               nq = as.integer(nY),
               xq = as.double(Q$x),
               yq = as.double(Q$y),
               nv = as.integer(Lvert$n),
               xv = as.double(Lvert$x),
               yv = as.double(Lvert$y),
               ns = as.double(L$n),
               from = as.integer(from0),
               to = as.integer(to0),
               dpath = as.double(dpath),
               psegmap = as.integer(Xsegmap),
               qsegmap = as.integer(Ysegmap),
               idP = as.integer(iX),
               idQ = as.integer(iY),
               huge = as.double(huge),
               nndist = as.double(nnd),
               nnwhich = as.integer(nnw))
    }
    nnd <- zz$nndist
    # convert C indexing to R indexing
    nnw <- zz$nnwhich + 1L
    # any zeroes occur if points have no neighbours.
    nnw[nnw == 0] <- NA
  }
  result <- data.frame(dist=nnd, which=nnw)[, what]
  return(result)
}
