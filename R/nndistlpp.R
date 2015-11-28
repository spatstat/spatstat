#
# nndistlpp.R
#
#  $Revision: 1.7 $ $Date: 2015/10/21 09:06:57 $
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
  L <- as.linnet(X$domain, sparse=FALSE)
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
    nseg <- length(from0)
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
             ns = as.integer(nseg),
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
  pro <- loco$seg %orifnull% nearestsegment(X, L$lines)

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
    nseg <- length(from0)
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
             ns = as.integer(nseg),
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

nncross.lpp <- function(X, Y, iX=NULL, iY=NULL,
                        what = c("dist", "which"), ..., method="C") {
  stopifnot(inherits(X, "lpp"))
  stopifnot(inherits(Y, "lpp"))
  what   <- match.arg(what, choices=c("dist", "which"), several.ok=TRUE)
  stopifnot(method %in% c("C", "interpreted"))
  if(is.null(iX) != is.null(iY))
    stop("If one of iX, iY is given, then both must be given")
  exclude <- (!is.null(iX) || !is.null(iY))

  check <- resolve.defaults(list(...), list(check=TRUE))$check
  if(check && !identical(as.linnet(X, sparse=TRUE),
                         as.linnet(Y, sparse=TRUE)))
    stop("X and Y are on different linear networks")

  fast <- (method == "C") && !exclude && spatstat.options("Cnncrosslpp")

  if(!fast) {
    ## require dpath matrix
    Xsparse <- identical(domain(X)$sparse, TRUE)
    Ysparse <- identical(domain(Y)$sparse, TRUE)
    L <- if(!Xsparse && Ysparse) as.linnet(X) else
         if(Xsparse && !Ysparse) as.linnet(Y) else
         as.linnet(X, sparse=FALSE)
  } else L <- as.linnet(X)
  #
  nX <- npoints(X)
  nY <- npoints(Y)
  P <- as.ppp(X)
  Q <- as.ppp(Y)
  #
  Lvert <- L$vertices
  from  <- L$from
  to    <- L$to
  if(fast) {
    seglengths <- lengths.psp(as.psp(L))
  } else {
    dpath <- L$dpath
  }
  
  # deal with null cases
  if(nX == 0)
    return(data.frame(dist=numeric(0), which=integer(0))[, what])
  if(nY == 0)
    return(data.frame(dist=rep(Inf, nX), which=rep(NA_integer_, nX))[, what])

  # find nearest segment for each point
  Xcoords <- coords(X)
  Ycoords <- coords(Y)
  Xpro <- Xcoords$seg
  Ypro <- Ycoords$seg

  # handle serial numbers
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
    nseg <- length(from0)
    Xsegmap <- Xpro - 1
    Ysegmap <- Ypro - 1
    # upper bound on interpoint distance
    huge <- if(!fast) {
      max(dpath) + 2 * diameter(Frame(L))
    } else {
      sum(seglengths)
    }
    # space for result
    nnd <- double(nX)
    nnw <- integer(nX)
    # call C
    if(fast) {
      ## experimental faster code
      message("Using experimental fast code for nncross.lpp")
      ooX <- order(Xsegmap)
      ooY <- order(Ysegmap)
      tol <- max(.Machine$double.eps,
                 diameter(Frame(L))/2^20)
      zz <- .C("linSnndwhich",
               np = as.integer(nX),
               sp = as.integer(Xsegmap[ooX]),
               tp = as.double(Xcoords$tp[ooX]),
               nq = as.integer(nY),
               sq = as.integer(Ysegmap[ooY]),
               tq = as.double(Ycoords$tp[ooY]),
               nv = as.integer(Lvert$n),
               ns = as.integer(nseg),
               from = as.integer(from0),
               to = as.integer(to0),
               seglen = as.double(seglengths), 
               huge = as.double(huge),
               tol = as.double(tol), 
               nndist = as.double(nnd),
               nnwhich = as.integer(nnw))
      zznd <- zz$nndist
      zznw <- zz$nnwhich + 1L
      if(any(notfound <- (zznw == 0))) {
        zznd[notfound] <- NA
        zznw[notfound] <- NA
      }
      nnd[ooX] <- zznd
      nnw[ooX] <- ooY[zznw]
    } else if(!exclude) {
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
               ns = as.integer(nseg),
               from = as.integer(from0),
               to = as.integer(to0),
               dpath = as.double(dpath),
               psegmap = as.integer(Xsegmap),
               qsegmap = as.integer(Ysegmap),
               huge = as.double(huge),
               nndist = as.double(nnd),
               nnwhich = as.integer(nnw))
      nnd <- zz$nndist
      nnw <- zz$nnwhich + 1L
    } else {
      ## excluding certain pairs
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
               ns = as.integer(nseg),
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
      nnd <- zz$nndist
      nnw <- zz$nnwhich + 1L
    }
    # any zeroes occur if points have no neighbours.
    nnw[nnw == 0] <- NA
  }
  result <- data.frame(dist=nnd, which=nnw)[, what]
  return(result)
}
