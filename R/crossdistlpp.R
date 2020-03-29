#
# crossdistlpp.R
#
#  $Revision: 1.11 $ $Date: 2020/03/29 09:08:10 $
#
#  crossdist.lpp
#        Calculates the shortest-path distance from each point of X
#        to each point of Y, where X and Y are point patterns
#        on the same linear network.
#

crossdist.lpp <- function(X, Y, ..., method="C", check=TRUE) {
  stopifnot(inherits(X, "lpp"))
  method <- match.arg(method, c("C", "interpreted"))
  L <- as.linnet(X)
  if(check) {
    LY <- as.linnet(Y) 
    if(!identical(L, LY))
      stop("X and Y are on different linear networks")
  }
  nX <- npoints(X)
  nY <- npoints(Y)

  crossdistmat <- matrix(Inf,nX,nY)
  
  if(!is.connected(L)) {
    #' disconnected network
    lab <- connected(L, what="labels")
    subsets <- split(seq_len(nvertices(L)), lab)
    for(subi in subsets) {
      Xi <- thinNetwork(X, retainvertices=subi)
      Yi <- thinNetwork(Y, retainvertices=subi)
      whichX <- attr(Xi, "retainpoints")      
      whichY <- attr(Yi, "retainpoints")      
      crossdistmat[whichX, whichY] <- crossdist.lpp(Xi, Yi, method=method)
    }
    return(crossdistmat)
  }

  ## ----------- network is connected ------------------------
  ## Extract network data
  Lvert <- L$vertices
  from  <- L$from
  to    <- L$to
  dpath <- L$dpath
  nseg  <- length(from)
  sparse <- L$sparse || is.null(dpath)
  ## Extract point coordinates
  P <- as.ppp(X)
  Q <- as.ppp(Y)
  ## local coordinates
  cooX <- coords(X, local=TRUE, spatial=FALSE, temporal=FALSE)
  cooY <- coords(Y, local=TRUE, spatial=FALSE, temporal=FALSE)
  Xseg <- cooX$seg
  Yseg <- cooY$seg
  ## 
  if(sparse) {
    ## new C code for sparse representation
    tX <- cooX$tp
    tY <- cooY$tp
    ## convert indices to start at 0
    from0 <- from - 1L
    to0   <- to - 1L
    Xseg0 <- Xseg - 1L
    Yseg0 <- Yseg - 1L
    ## sort each set of points by increasing segment index
    ordX <- order(Xseg0, tX)
    Xseg0 <- Xseg0[ordX]
    tX    <- tX[ordX]
    ordY <- order(Yseg0, tY)
    Yseg0 <- Yseg0[ordY]
    tY    <- tY[ordY]
    ## network info
    seglen <- lengths_psp(L$lines)
    huge <- 2 * sum(seglen)
    tol <- L$toler %orifnull% default.linnet.tolerance(L)
    ## 
    zz <- .C("linScrossdist",
             np = as.integer(nX),
             sp = as.integer(Xseg0),
             tp = as.double(tX),
             nq = as.integer(nY),
             sq = as.integer(Yseg0),
             tq = as.double(tY),
             nv = as.integer(Lvert$n),
             ns = as.integer(nseg),
             from   = as.integer(from0),
             to     = as.integer(to0),
             seglen = as.double(seglen),
             huge   = as.double(huge),
             tol    = as.double(tol),
             dist   = as.double(numeric(nX * nY)),
             PACKAGE = "spatstat")
    crossdistmat[ordX, ordY] <- zz$dist
  } else {
    switch(method,
           C = {
             ## older C code requiring non-sparse representation
             ## convert indices to start at 0
             from0 <- from - 1L
             to0   <- to - 1L
             Xsegmap <- Xseg - 1L
             Ysegmap <- Yseg - 1L
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
                      ns = as.integer(nseg),
                      from = as.integer(from0),
                      to = as.integer(to0),
                      dpath = as.double(dpath),
                      psegmap = as.integer(Xsegmap),
                      qsegmap = as.integer(Ysegmap),
                      answer = as.double(numeric(nX * nY)),
                      PACKAGE = "spatstat")
             crossdistmat <- matrix(zz$answer, nX, nY)
           },
           interpreted = {
             #' interpreted code
             #' loop through all pairs of data points
             for (i in 1:nX) {
               Xsegi <- Xseg[i]
               Xi <- P[i]
               nbi1 <- from[Xsegi]
               nbi2 <- to[Xsegi]
               vi1 <- Lvert[nbi1]
               vi2 <- Lvert[nbi2]   
               dXi1 <- crossdist(Xi, vi1)
               dXi2 <- crossdist(Xi, vi2)
               for (j in 1:nY) {
                 Yj <- Q[j]
                 Ysegj <- Yseg[j]
                 if(Xsegi == Ysegj) {
                   #' points i and j lie on the same segment
                   #' use Euclidean distance
                   d <- crossdist(Xi, Yj)
                 } else {
                   #' shortest path from i to j passes through ends of segments
                   nbj1 <- from[Ysegj]
                   nbj2 <- to[Ysegj]
                   vj1 <- Lvert[nbj1]
                   vj2 <- Lvert[nbj2]
                   #' Calculate shortest of 4 possible paths from i to j
                   d1Yj <- crossdist(vj1,Yj)
                   d2Yj <- crossdist(vj2,Yj)
                   d11 <- dXi1 + dpath[nbi1,nbj1] + d1Yj
                   d12 <- dXi1 + dpath[nbi1,nbj2] + d2Yj
                   d21 <- dXi2 + dpath[nbi2,nbj1] + d1Yj
                   d22 <- dXi2 + dpath[nbi2,nbj2] + d2Yj
                   d <- min(d11,d12,d21,d22)
                 }
                 #' store result
                 crossdistmat[i,j] <- d
               }
             }
           }
      )
  }
  return(crossdistmat)
}
