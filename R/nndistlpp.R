#
# nndistlpp.R
#
#  $Revision: 1.27 $ $Date: 2020/04/24 03:52:39 $
#
# Methods for nndist, nnwhich, nncross for linear networks
#
# nndist.lpp
#   Calculates the nearest neighbour distances in the shortest-path metric
#   for a point pattern on a linear network.

nndist.lpp <- function(X, ..., k=1, by=NULL, method="C") {
  stopifnot(inherits(X, "lpp"))
  stopifnot(method %in% c("C", "interpreted"))
  n <- npoints(X)
  k <- as.integer(k)
  stopifnot(all(k > 0))
  kmax <- max(k)

  if(!is.null(by)) 
    return(genericNNdistBy(X, by, k=k))

  L <- as.linnet(X)
  if(is.null(br <- L$boundingradius) || is.infinite(br)) {
    # network may be disconnected
    lab <- connected(L, what="labels")
    if(length(levels(lab)) > 1L) {
      # network is disconnected
      result <- matrix(Inf, n, length(k))
      # handle each connected component separately
      subsets <- split(seq_len(nvertices(L)), lab)
      for(i in seq_along(subsets)) {
        Xi <- thinNetwork(X, retainvertices=subsets[[i]])
        relevant <- attr(Xi, "retainpoints")      
        result[relevant, ] <- nndist.lpp(Xi, k=k, method=method)
      }
      return(result)
    }
  }
  
  toomany <- (kmax >= n-1)
  if(toomany) {
    ## not enough points to define kmax nearest neighbours
    result <- matrix(Inf, nrow=n, ncol=kmax)
    if(n <= 1) return(result[,k,drop=TRUE])
    ## reduce kmax to feasible value
    kmax <- n-1
    kuse <- k[k <= kmax]
  } else {
    kuse <- k
  }
  
  Y <- as.ppp(X)
  sparse <- identical(L$sparse, TRUE)

  ## find nearest segment for each point
  ## This is given by local coordinates, if available (spatstat >= 1.28-0)
  loco <- coords(X, local=TRUE, spatial=FALSE, temporal=FALSE)
  pro <- if(!is.null(seg <- loco$seg)) seg else nearestsegment(X, Lseg)
  
  if(method == "interpreted") {
    ## interpreted code 
    D <- pairdist(X, method="interpreted")
    diag(D) <- Inf
    ans <- if(kmax == 1) apply(D, 1, min) else
           t(apply(D, 1, orderstats, k=kuse))[,,drop=TRUE]
  } else if(!sparse && kmax == 1) {
    # C code for non-sparse network
    Lseg  <- L$lines
    Lvert <- L$vertices
    from  <- L$from
    to    <- L$to
    dpath <- L$dpath
    # convert indices to start at 0
    from0 <- from - 1L
    to0   <- to - 1L
    segmap <- pro - 1L
    nseg <- length(from0)
    # upper bound on interpoint distance
    huge <- max(dpath) + 2 * max(lengths_psp(Lseg))
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
             answer = as.double(ans),
             PACKAGE = "spatstat")
    ans <- zz$answer
  } else if(spatstat.options('Cnndistlpp')) {
    ## use new C routine
    Lseg  <- L$lines
    Lvert <- L$vertices
    from  <- L$from
    to    <- L$to
    ##
    nseg <- length(from)
    seglen <- lengths_psp(Lseg)
    ## convert indices to start at 0
    from0 <- from - 1L
    to0   <- to - 1L
    segmap <- pro - 1L
    tp <- loco$tp
    ## sort by segment index
    oo <- order(segmap, tp)
    segmap <- segmap[oo]
    tp <- tp[oo]
    # upper bound on interpoint distance
    huge <- sum(seglen)
    #' numerical tolerance
    tol <- max(.Machine$double.eps,
               diameter(Frame(L))/2^20)
    #'
    kmax1 <- kmax + 1L
    zz <- .C("linknnd",
             kmax = as.integer(kmax1),
             np = as.integer(n),
             sp = as.integer(segmap),
             tp = as.double(tp), 
             nv = as.integer(Lvert$n),
             ns = as.integer(nseg),
             from = as.integer(from0),
             to = as.integer(to0),
             seglen = as.double(seglen),
             huge = as.double(huge),
             tol = as.double(tol),
             nndist = as.double(numeric(n * kmax1)),
             nnwhich = as.integer(integer(n * kmax1)),
             PACKAGE = "spatstat")
    ans <- matrix(, n, kmax1)
    ans[oo, ] <- matrix(zz$nndist, n, kmax1, byrow=TRUE)
    # drop first column which is zero corresponding to j = i
    ans <- ans[, -1, drop=FALSE]
    colnames(ans) <- paste0("dist.", 1:ncol(ans))
    ans <- ans[,kuse]
  } else {    
    ## use fast code for nncross
    ans <- nncross(X, X, what="dist", k=kuse+1)
    if(is.matrix(ans) || is.data.frame(ans))
      colnames(ans) <- paste0("dist.", kuse)
  }
  if(!is.null(dim(ans))) {
    ans <- as.matrix(ans)
    rownames(ans) <- NULL
  } 
  if(!toomany)
    return(ans)
  result[, kuse] <- as.matrix(ans)
  colnames(result) <- paste0("dist.", 1:ncol(result))
  return(result[,k])
}

# nnwhich.lpp
# Identifies the nearest neighbours in the shortest-path metric
# for a point pattern on a linear network.
#

nnwhich.lpp <- function(X, ..., k=1, method="C") {
  stopifnot(inherits(X, "lpp"))
  stopifnot(method %in% c("C", "interpreted"))

  k <- as.integer(k)
  stopifnot(all(k > 0))
  kmax <- max(k)

  n <- npoints(X)

  L <- as.linnet(X)
  if(is.null(br <- L$boundingradius) || is.infinite(br)) {
    # network may be disconnected
    lab <- connected(L, what="labels")
    if(length(levels(lab)) > 1L) {
      # network is disconnected
      result <- matrix(NA_integer_, n, length(k))
      # handle each connected component separately
      subsets <- split(seq_len(nvertices(L)), lab)
      for(i in seq_along(subsets)) {
        Xi <- thinNetwork(X, retainvertices=subsets[[i]])
        relevant <- attr(Xi, "retainpoints")      
        result[relevant, ] <- nnwhich.lpp(Xi, k=k, method=method)
      }
      return(result)
    }
  }
  
  toomany <- (kmax >= n-1)
  if(toomany) {
    ## not enough points to define kmax nearest neighbours
    result <- matrix(NA_integer_, nrow=n, ncol=kmax)
    if(n <= 1) return(result[,k,drop=TRUE])
    ## reduce kmax to feasible value
    kmax <- n-1
    kuse <- k[k <= kmax]
  } else {
    kuse <- k
  }
  
  #
  Y <- as.ppp(X)
  sparse <- identical(L$sparse, TRUE)
  
  ## find nearest segment for each point
  ## This is given by local coordinates, if available (spatstat >= 1.28-0)
  loco <- coords(X, local=TRUE, spatial=FALSE, temporal=FALSE)
  pro <- if(!is.null(seg <- loco$seg)) seg else nearestsegment(X, Lseg)

  if(method == "interpreted") {
    D <- pairdist(X, method="interpreted")
    diag(D) <- Inf
    nnw <- if(kmax == 1) apply(D, 1, which.min) else
           t(apply(D, 1, orderwhich, k=kuse))[,,drop=TRUE]
  } else if(!sparse && kmax == 1) {
    # C code for non-sparse network
    ##
    Lseg  <- L$lines
    Lvert <- L$vertices
    from  <- L$from
    to    <- L$to
    dpath <- L$dpath
    ## convert indices to start at 0
    from0 <- from - 1L
    to0   <- to - 1L
    segmap <- pro - 1L
    nseg <- length(from0)
    # upper bound on interpoint distance
    huge <- max(dpath) + 2 * max(lengths_psp(Lseg))
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
             nnwhich = as.integer(nnw),
             PACKAGE = "spatstat")
    # convert C indexing to R indexing
    nnw <- zz$nnwhich + 1L
    # any zeroes occur if points have no neighbours.
    nnw[nnw == 0] <- NA
  } else if(spatstat.options('Cnndistlpp')) {
    ## use new C routine
    Lseg  <- L$lines
    Lvert <- L$vertices
    from  <- L$from
    to    <- L$to
    ##
    nseg <- length(from)
    seglen <- lengths_psp(Lseg)
    ## convert indices to start at 0
    from0 <- from - 1L
    to0   <- to - 1L
    segmap <- pro - 1L
    tp <- loco$tp
    ## sort by segment index
    oo <- order(segmap, tp)
    segmap <- segmap[oo]
    tp <- tp[oo]
    # upper bound on interpoint distance
    huge <- sum(seglen)
    #' numerical tolerance
    tol <- max(.Machine$double.eps,
               diameter(Frame(L))/2^20)
    #'
    kmax1 <- kmax + 1L
    zz <- .C("linknnd",
             kmax = as.integer(kmax1),
             np = as.integer(n),
             sp = as.integer(segmap),
             tp = as.double(tp), 
             nv = as.integer(Lvert$n),
             ns = as.integer(nseg),
             from = as.integer(from0),
             to = as.integer(to0),
             seglen = as.double(seglen),
             huge = as.double(huge),
             tol = as.double(tol),
             nndist = as.double(numeric(n * kmax1)),
             nnwhich = as.integer(integer(n * kmax1)),
             PACKAGE = "spatstat")
    nnw <- matrix(, n, kmax1)
    nnw[oo, ] <- matrix(oo[zz$nnwhich + 1L], n, kmax1, byrow=TRUE)
    # drop first column which is j = i
    nnw <- nnw[, -1, drop=FALSE]
    colnames(nnw) <- paste0("which.", 1:ncol(nnw))
    nnw <- nnw[,kuse]
  } else {
    ## use fast code for nncross
    nnw <- nncross(X, X, what="which", k=kuse+1)
    if(is.matrix(nnw) || is.data.frame(nnw))
      colnames(nnw) <- paste0("which.", kuse)
  }
  if(!is.null(dim(nnw))) {
    nnw <- as.matrix(nnw)
    rownames(nnw) <- NULL
  }
  if(!toomany)
    return(nnw)
  result[, kuse] <- as.matrix(nnw)
  colnames(result) <- paste0("which.", 1:ncol(result))
  return(result[,k])
}

# nncross.lpp
# Identifies the nearest neighbours in the shortest-path metric
# from one point pattern on a linear network to ANOTHER pattern
# on the SAME network.
#

nncross.lpp <- local({

  nncross.lpp <- function(X, Y, iX=NULL, iY=NULL,
                          what = c("dist", "which"), ..., k=1, method="C") {
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

  # internal use only
  format <- resolve.defaults(list(...), list(format="data.frame"))$format

  nX <- npoints(X)
  nY <- npoints(Y)

  L <- domain(X)
  if(is.null(br <- L$boundingradius) || is.infinite(br)) {
    # network may be disconnected
    lab <- connected(L, what="labels")
    if(length(levels(lab)) > 1L) {
      # network is disconnected
      # handle each connected component separately
      subsets <- split(seq_len(nvertices(L)), lab)
      nndistmat <- if("dist" %in% what) matrix(Inf, nX, length(k)) else NULL
      nnwhichmat <-
         if("which" %in% what) matrix(NA_integer_, nX, length(k)) else NULL
      for(i in seq_along(subsets)) {
        subi <- subsets[[i]]
        Xi <- thinNetwork(X, retainvertices=subi)
        useX <- attr(Xi, "retainpoints")      
        Yi <- thinNetwork(Y, retainvertices=subi)
        useY <- attr(Yi, "retainpoints")
	z <- nncross.lpp(Xi, Yi,
	                 iX = iX[useX], iY=iY[useY],
	                 what=what, k=k, method=method,
			 format="list")
        if("dist" %in% what)
	   nndistmat[useX, ] <- z$dist
        if("which" %in% what)
	   nnwhichmat[useX, ] <- which(useY)[z$which]
      }
      return(shapedresult(dist=nndistmat,
                          which=nnwhichmat,
                          what=what, format=format))
    }
  }

  koriginal <- k <- as.integer(k)
  stopifnot(all(k > 0))
  kmax <- max(k)

  #' decide which algorithm to use
  #' fast C algorithm 
  fast <- (method == "C") && (spatstat.options("Cnncrosslpp") || (kmax > 1))
  #' slower C algorithm for exclusion case for k=1
  excludeinC <- exclude && (method == "C") && !fast && (k == 1)
  excludeinR <- exclude && !excludeinC

  if(excludeinR) {
    #' compute k+1 neighbours in C, then filter in R
    kmax <- kmax+1
    k <- 1:kmax
  }
  
  toomany <- (kmax > nY)
  if(toomany) {
    paddist <- matrix(Inf, nX, kmax)
    padwhich <- matrix(NA_integer_, nX, kmax)
    kmax <- nY
    kuse <- k[k <= kmax]
  } else {
    kuse <- k
  }

  if(length(kuse) == 0) {
    # None of the required values are defined
    nnd <- paddist
    nnw <- padwhich
    maxk <- max(k)
    colnames(nnd) <- paste0("dist.", seq_len(maxk))
    colnames(nnd) <- paste0("dist.", seq_len(maxk))
    nnd <- nnd[,k,drop=TRUE]
    nnw <- nnw[,k,drop=TRUE]
    return(shapedresult(dist=nnd, which=nnw, what=what, format=format))
  }
  
  need.dist <- ("dist" %in% what) || excludeinR
  need.which <- ("which" %in% what) || excludeinR
  
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
    seglengths <- lengths_psp(as.psp(L))
  } else {
    dpath <- L$dpath
  }
  
  # deal with null cases
  if(nX == 0) 
    return(shapedresult(dist=numeric(0),
                        which=integer(0),
                        what=what, format=format))

  if(nY == 0)
    return(shapedresult(dist=rep(Inf, nX),
                        which=rep(NA_integer_, nX),
                        what=what, format=format))

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
    ## interpreted code
    D <- crossdist(X, Y, method="interpreted")
    if(exclude)
      D[outer(iX, iY, "==")] <- Inf
    nnd <- nnw <- NULL
    if(need.dist) {
      nnd <- if(kmax == 1) apply(D, 1, min) else
             t(apply(D, 1, orderstats, k=kuse))[,,drop=TRUE]
    }
    if(need.which) {
      nnw <- if(kmax == 1) apply(D, 1, which.min) else
             t(apply(D, 1, orderwhich, k=kuse))[,,drop=TRUE]
    } 
  } else {
    ## C code
    ## convert indices to start at 0
    from0 <- from - 1L
    to0   <- to - 1L
    nseg <- length(from0)
    Xsegmap <- Xpro - 1L
    Ysegmap <- Ypro - 1L
    ## upper bound on interpoint distance
    huge <- if(!fast) {
      max(dpath) + 2 * diameter(Frame(L))
    } else {
      sum(seglengths)
    }
    ## space for result
    nnd <- double(nX * kmax)
    nnw <- integer(nX * kmax)
    ## call C
    if(fast) {
      ## experimental faster code
      ooX <- order(Xsegmap)
      ooY <- order(Ysegmap)
      tol <- max(.Machine$double.eps,
                 diameter(Frame(L))/2^20)
      if(kmax > 1) {
        zz <- .C("linknncross",
                 kmax = as.integer(kmax),
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
                 nnwhich = as.integer(nnw),
                 PACKAGE = "spatstat")
        zznd <- matrix(zz$nndist, ncol=kmax, byrow=TRUE)
        zznw <- matrix(zz$nnwhich + 1L, ncol=kmax, byrow=TRUE)
        if(any(notfound <- (zznw == 0))) {
          zznd[notfound] <- NA
          zznw[notfound] <- NA
        }
        nnd <- matrix(nnd, nX, kmax)
        nnw <- matrix(nnw, nX, kmax)
        nnd[ooX, ] <- zznd
        nnw[ooX, ] <- ooY[zznw]
        colnames(nnd) <- colnames(nnw) <- seq_len(kmax)
        if(!identical(kuse, seq_len(kmax))) {
          nnd <- nnd[,kuse,drop=FALSE]
          nnw <- nnw[,kuse,drop=FALSE]
          if(length(kuse) == 1) {
            colnames(nnd) <- paste0("dist.", kuse)
            colnames(nnw) <- paste0("which.", kuse)
          }
        }
      } else {
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
                 nnwhich = as.integer(nnw),
                 PACKAGE = "spatstat")
        zznd <- zz$nndist
        zznw <- zz$nnwhich + 1L
        if(any(notfound <- (zznw == 0))) {
          zznd[notfound] <- NA
          zznw[notfound] <- NA
        }
        nnd[ooX] <- zznd
        nnw[ooX] <- ooY[zznw]
      }
    } else {
      ## slower code requiring dpath matrix
      if(!excludeinC) {
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
                 nnwhich = as.integer(nnw),
                 PACKAGE = "spatstat")
        nnd <- zz$nndist
        nnw <- zz$nnwhich + 1L
      } else {
        ## excluding certain pairs (k=1)
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
                 nnwhich = as.integer(nnw),
                 PACKAGE = "spatstat")
        nnd <- zz$nndist
        nnw <- zz$nnwhich + 1L
      }
      ## any zeroes occur if points have no neighbours.
      nnw[nnw == 0] <- NA
    }
  }
  
  if(toomany) {
    ## Nearest neighbours were undefined for some large values of k.
    ## Insert results obtained for valid 'k' back into matrix of NA/Inf
    if(need.dist) {
      paddist[,kuse] <- as.matrix(nnd)
      nnd <- paddist
    }
    if(need.which) {
      padwhich[,kuse] <- as.matrix(nnw)
      nnw <- padwhich
    }
  }
  if(excludeinR) {
    ## now find neighbours that don't have the same id number
    if(!is.matrix(nnw)) nnw <- as.matrix(nnw, ncol=1)
    if(!is.matrix(nnd)) nnd <- as.matrix(nnd, ncol=1)
    avoid <- matrix(iX[as.vector(row(nnw))] != iY[as.vector(nnw)],
                    nrow=nrow(nnw), ncol=ncol(nnw))
    colind <- apply(avoid, 1, whichcoltrue, m=seq_len(ncol(avoid)-1))
    colind <- if(is.matrix(colind)) t(colind) else matrix(colind, ncol=1)
    rowcol <- cbind(as.vector(row(colind)), as.vector(colind))
    nnd <- matrix(nnd[rowcol], nrow=nX)
    nnw <- matrix(nnw[rowcol], nrow=nX)
    nnd <- nnd[,koriginal]
    nnw <- nnw[,koriginal]
  }
  return(shapedresult(dist=nnd, which=nnw, what=what, format=format))
}

  whichcoltrue <- function(x, m) which(x)[m]

  shapedresult <- function(dist, which,
                           what=c("dist", "which"),
                           format="data.frame") {
    #' idiom to return result in correct format
    result <- list(dist=dist, which=which)[what]
    if(format == "data.frame")
      result <- as.data.frame(result)[,,drop=TRUE]
    return(result)
  }
  
  nncross.lpp
})
