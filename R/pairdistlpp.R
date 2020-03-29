#
# pairdistlpp.R
#
#  $Revision: 1.20 $ $Date: 2020/03/29 09:08:05 $
#
#
#  pairdist.lpp
#        Calculates the shortest-path distance between each pair of points
#        in a point pattern on a linear network.
#

pairdist.lpp <- function(X, ..., method="C") {
  stopifnot(inherits(X, "lpp"))
  method <- match.arg(method, c("C", "interpreted", "testsymm"))
  #
  n <- npoints(X)
  pairdistmat <- matrix(Inf,n,n)
  diag(pairdistmat) <- 0
  if(n <= 1) return(pairdistmat)
  #
  L <- as.linnet(X)
  #
  if(!is.connected(L)) {
    #' disconnected network
    lab <- connected(L, what="labels")
    subsets <- split(seq_len(nvertices(L)), lab)
    for(i in seq_along(subsets)) {
      Xi <- thinNetwork(X, retainvertices=subsets[[i]])
      witch <- attr(Xi, "retainpoints")      
      pairdistmat[witch, witch] <- pairdist.lpp(Xi, method=method)
    }
    return(pairdistmat)
  }
  # 
  ## ----------- network is connected ------------------------
  ## Extract network data
  Lvert <- L$vertices
  nvert <- npoints(Lvert)
  from  <- L$from
  to    <- L$to
  dpath <- L$dpath
  nseg  <- length(from)
  sparse <- L$sparse || is.null(dpath)
  ## Point coordinates
  Y <- as.ppp(X)
  cooX <- coords(X, local=TRUE, spatial=FALSE, temporal=FALSE)
  Xseg <- cooX$seg
  ##
  if(sparse || method == "testsymm") {
    ## new C code for sparse representation
    tX <- cooX$tp
    ## convert indices to start at 0
    from0 <- from - 1L
    to0   <- to - 1L
    Xseg0 <- Xseg - 1L
    ## sort points by increasing segment index
    ordX <- order(Xseg0, tX)
    Xseg0 <- Xseg0[ordX]
    tX    <- tX[ordX]
    ## network info
    seglen <- lengths_psp(L$lines)
    huge <- 2 * sum(seglen)
    tol <- L$toler %orifnull% default.linnet.tolerance(L)
    ##
    if(method == "testsymm") {
      ## test whether the sparse algorithm yields symmetric results
      zz <- .C("linSpairUdist",
               np = as.integer(n),
               sp = as.integer(Xseg0),
               tp = as.double(tX),
               nv = as.integer(nvert),
               ns = as.integer(nseg),
               from   = as.integer(from0),
               to     = as.integer(to0),
               seglen = as.double(seglen),
               huge   = as.double(huge),
               tol    = as.double(tol),
               dist   = as.double(numeric(n * n)),
               PACKAGE = "spatstat")
    } else {
      ## use algorithm which exploits symmetry
      zz <- .C("linSpairdist",
               np = as.integer(n),
               sp = as.integer(Xseg0),
               tp = as.double(tX),
               nv = as.integer(nvert),
               ns = as.integer(nseg),
               from   = as.integer(from0),
               to     = as.integer(to0),
               seglen = as.double(seglen),
               huge   = as.double(huge),
               tol    = as.double(tol),
               dist   = as.double(numeric(n * n)),
               PACKAGE = "spatstat")
    }
    pairdistmat[ordX, ordX] <- zz$dist
  } else {
    switch(method,
           interpreted = {
             ## loop through all pairs of data points
             for (i in 1:(n-1)) {
               Xsegi <- Xseg[i]
               Xi <- Y[i]
               nbi1 <- from[Xsegi]
               nbi2 <- to[Xsegi]
               vi1 <- Lvert[nbi1]
               vi2 <- Lvert[nbi2]   
               dXi1 <- crossdist(Xi, vi1)
               dXi2 <- crossdist(Xi, vi2)
               for (j in (i+1):n) {
                 Xj <- Y[j]
                 Xsegj <- Xseg[j]
                 if(Xsegi == Xsegj) {
                   ## points i and j lie on the same segment
                   ## use Euclidean distance
                   d <- crossdist(Xi, Xj)
                 } else {
                   ## shortest path from i to j passes through ends of segments
                   nbj1 <- from[Xsegj]
                   nbj2 <- to[Xsegj]
                   vj1 <- Lvert[nbj1]
                   vj2 <- Lvert[nbj2]
                   ## Calculate shortest of 4 possible paths from i to j
                   d1Xj <- crossdist(vj1,Xj)
                   d2Xj <- crossdist(vj2,Xj)
                   d11 <- dXi1 + dpath[nbi1,nbj1] + d1Xj
                   d12 <- dXi1 + dpath[nbi1,nbj2] + d2Xj
                   d21 <- dXi2 + dpath[nbi2,nbj1] + d1Xj
                   d22 <- dXi2 + dpath[nbi2,nbj2] + d2Xj
                   d <- min(d11,d12,d21,d22)
                 }
                 ## store result
                 pairdistmat[i,j] <- pairdistmat[j,i] <- d
               }
             }
           },
           C = {
             ## C code using non-sparse representation
             ## convert indices to start at 0
             from0 <- from - 1L
             to0   <- to - 1L
             segmap <- Xseg - 1L
             zz <- .C("linpairdist",
                      np = as.integer(n),
                      xp = as.double(Y$x),
                      yp = as.double(Y$y),
                      nv = as.integer(nvert),
                      xv = as.double(Lvert$x),
                      yv = as.double(Lvert$y),
                      ns = as.double(L$n),
                      from = as.integer(from0),
                      to = as.integer(to0),
                      dpath = as.double(dpath),
                      segmap = as.integer(segmap),
                      answer = as.double(numeric(n*n)),
                      PACKAGE = "spatstat")
             pairdistmat <- matrix(zz$answer, n, n)
           })
  }
  return(pairdistmat)
}

