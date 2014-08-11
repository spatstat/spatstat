#
# pairdistlpp.R
#
#  $Revision: 1.9 $ $Date: 2012/10/13 02:25:43 $
#
#
#  pairdist.lpp
#        Calculates the shortest-path distance between each pair of points
#        in a point pattern on a linear network.
#

pairdist.lpp <- function(X, ..., method="C") {
  stopifnot(inherits(X, "lpp"))
  stopifnot(method %in% c("C", "interpreted"))
  #
  L <- as.linnet(X)
  Y <- as.ppp(X)
  n <- npoints(Y)
  #
  Lseg  <- L$lines
  Lvert <- L$vertices
  from  <- L$from
  to    <- L$to
  dpath <- L$dpath
  
  # nearest segment for each point
  pro <- coords(X, local=TRUE, spatial=FALSE, temporal=FALSE)$seg

  pairdistmat <- matrix(0,n,n)

  if(method == "interpreted") {
    # loop through all pairs of data points
    for (i in 1:(n-1)) {
      proi <- pro[i]
      Xi <- Y[i]
      nbi1 <- from[proi]
      nbi2 <- to[proi]
      vi1 <- Lvert[nbi1]
      vi2 <- Lvert[nbi2]   
      dXi1 <- crossdist(Xi, vi1)
      dXi2 <- crossdist(Xi, vi2)
      for (j in (i+1):n) {
        Xj <- Y[j]
        proj <- pro[j]
        if(proi == proj) {
          # points i and j lie on the same segment
          # use Euclidean distance
          d <- crossdist(Xi, Xj)
        } else {
          # shortest path from i to j passes through ends of segments
          nbj1 <- from[proj]
          nbj2 <- to[proj]
          vj1 <- Lvert[nbj1]
          vj2 <- Lvert[nbj2]
          # Calculate shortest of 4 possible paths from i to j
          d1Xj <- crossdist(vj1,Xj)
          d2Xj <- crossdist(vj2,Xj)
          d11 <- dXi1 + dpath[nbi1,nbj1] + d1Xj
          d12 <- dXi1 + dpath[nbi1,nbj2] + d2Xj
          d21 <- dXi2 + dpath[nbi2,nbj1] + d1Xj
          d22 <- dXi2 + dpath[nbi2,nbj2] + d2Xj
          d <- min(d11,d12,d21,d22)
        }
        # store result
        pairdistmat[i,j] <- pairdistmat[j,i] <- d
      }
    }
  } else {
    # C code
    # convert indices to start at 0
    from0 <- from - 1L
    to0   <- to - 1L
    segmap <- pro - 1L
    zz <- .C("linpairdist",
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
             answer = as.double(pairdistmat),
             PACKAGE="spatstat")
    pairdistmat <- matrix(zz$answer, n, n)
  }
  return(pairdistmat)
}

