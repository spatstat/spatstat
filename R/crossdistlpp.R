#
# crossdistlpp.R
#
#  $Revision: 1.6 $ $Date: 2017/06/05 10:31:58 $
#
#  crossdist.lpp
#        Calculates the shortest-path distance from each point of X
#        to each point of Y, where X and Y are point patterns
#        on the same linear network.
#

crossdist.lpp <- function(X, Y, ..., method="C") {
  stopifnot(inherits(X, "lpp"))
  stopifnot(method %in% c("C", "interpreted"))
  check <- resolve.defaults(list(...), list(check=TRUE))$check
  #
  nX <- npoints(X)
  nY <- npoints(Y)
  #
  L <- as.linnet(X, sparse=FALSE)
  if(check) {
    LY <- as.linnet(Y, sparse=FALSE)
    if(!identical(L, LY))
      stop("X and Y are on different linear networks")
  }

  if(any(is.infinite(L$dpath))) {
    #' disconnected network
    lab <- connected(L, what="labels")
    subsets <- split(seq_len(nvertices(L)), lab)
    crossdistmat <- matrix(Inf,nX,nY)
    for(subi in subsets) {
      Xi <- thinNetwork(X, retainvertices=subi)
      Yi <- thinNetwork(Y, retainvertices=subi)
      whichX <- attr(Xi, "retainpoints")      
      whichY <- attr(Yi, "retainpoints")      
      crossdistmat[whichX, whichY] <- crossdist.lpp(Xi, Yi, method=method)
    }
    return(crossdistmat)
  }

  # network is connected
  
  P <- as.ppp(X)
  Q <- as.ppp(Y)
  #
#  Lseg  <- L$lines
  Lvert <- L$vertices
  from  <- L$from
  to    <- L$to
  dpath <- L$dpath
  
  # nearest segment for each point
  Xpro <- coords(X, local=TRUE, spatial=FALSE, temporal=FALSE)$seg
  Ypro <- coords(Y, local=TRUE, spatial=FALSE, temporal=FALSE)$seg

  if(method == "interpreted") {
    # loop through all pairs of data points
    crossdistmat <- matrix(,nX,nY)
    for (i in 1:nX) {
      Xproi <- Xpro[i]
      Xi <- P[i]
      nbi1 <- from[Xproi]
      nbi2 <- to[Xproi]
      vi1 <- Lvert[nbi1]
      vi2 <- Lvert[nbi2]   
      dXi1 <- crossdist(Xi, vi1)
      dXi2 <- crossdist(Xi, vi2)
      for (j in 1:nY) {
        Yj <- Q[j]
        Yproj <- Ypro[j]
        if(Xproi == Yproj) {
          # points i and j lie on the same segment
          # use Euclidean distance
          d <- crossdist(Xi, Yj)
        } else {
          # shortest path from i to j passes through ends of segments
          nbj1 <- from[Yproj]
          nbj2 <- to[Yproj]
          vj1 <- Lvert[nbj1]
          vj2 <- Lvert[nbj2]
          # Calculate shortest of 4 possible paths from i to j
          d1Yj <- crossdist(vj1,Yj)
          d2Yj <- crossdist(vj2,Yj)
          d11 <- dXi1 + dpath[nbi1,nbj1] + d1Yj
          d12 <- dXi1 + dpath[nbi1,nbj2] + d2Yj
          d21 <- dXi2 + dpath[nbi2,nbj1] + d1Yj
          d22 <- dXi2 + dpath[nbi2,nbj2] + d2Yj
          d <- min(d11,d12,d21,d22)
        }
        # store result
        crossdistmat[i,j] <- d
      }
    }
  } else {
    # C code
    # convert indices to start at 0
    from0 <- from - 1L
    to0   <- to - 1L
    Xsegmap <- Xpro - 1L
    Ysegmap <- Ypro - 1L
    zz <- .C("lincrossdist",
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
             answer = as.double(numeric(nX * nY)),
             PACKAGE = "spatstat")
    crossdistmat <- matrix(zz$answer, nX, nY)
  }
  return(crossdistmat)
}
