#'  nnfromvertex.R
#'
#'  Nearest data point to each vertex of a network
#'
#'  $Revision: 1.3 $  $Date: 2020/03/16 10:28:51 $
#'

nnfromvertex <- function(X, what=c("dist", "which"), k=1) {
  stopifnot(is.lpp(X))
  what <- match.arg(what, several.ok=TRUE)
  nX <- npoints(X)
  nv <- nvertices(domain(X))

  #' k can be a single integer or an integer vector
  if(length(k) == 0)
    stop("k is an empty vector")
  else if(length(k) == 1) {
    if(k != round(k) || k <= 0)
      stop("k is not a positive integer")
  } else {
    if(any(k != round(k)) || any(k <= 0))
      stop(paste("some entries of the vector",
           sQuote("k"), "are not positive integers"))
  }
  k <- as.integer(k)
  kmax <- max(k)

  #' Initialise results
  nnd <- matrix(Inf,         nrow=nv, ncol=kmax)
  nnw <- matrix(NA_integer_, nrow=nv, ncol=kmax)
  colnames(nnd) <- colnames(nnw) <- 1:kmax
  
  #' Trivial cases
  if(nX > 0) {
    #' Unique points, remembering original sequence
    ii <- which(!duplicated(X))
    uX <- X[ii]
    #' local coordinates
    coUX <- coords(uX)[, c("seg", "tp")]
    #' add label from original sequence index
    coUX$lab <- ii
    #' reorder
    oo <- with(coUX, order(seg, tp))
    coUXord <- coUX[oo, , drop=FALSE]
    seg <- coUXord$seg
    tp  <- coUXord$tp
    #' network data
    L <- domain(X)
    nv <- nvertices(L)
    ns <- nsegments(L)
    seglen <- lengths_psp(as.psp(L))
    from <- L$from
    to   <- L$to
    #' upper bound on interpoint distance
    huge <- sum(seglen)
    #' numerical tolerance for nnwhich
    tol <- max(sqrt(.Machine$double.eps), diameter(Frame(L))/2^20)
    #' ..............................................
    #' number of neighbours that are well-defined
    kmaxcalc <- min(nX, kmax)
    #' calculate k-nn distances and identifiers for 1 <= k <= kmaxcalc
    z <- vnnFind(seg, tp, ns, nv, from, to, seglen, huge, tol, kmax=kmaxcalc)
    vnndist  <- z$vnndist
    vnnwhich <- z$vnnwhich
    #' map identifiers back to original data pattern
    vnnwhich <- coUXord$lab[vnnwhich]
    #' insert results in correct places
    nnd[, 1:kmaxcalc] <- vnndist
    nnw[, 1:kmaxcalc] <- vnnwhich
  }
  #' extract required values
  nnd <- nnd[,k, drop=TRUE]
  nnw <- nnw[,k, drop=TRUE]
  if(identical(what, "dist")) return(nnd)
  if(identical(what, "which")) return(nnw)
  return(cbind(data.frame(dist=nnd), data.frame(which=nnw)))
}


vnnFind <- function(seg, tp, ns, nv, from, to, seglen, huge, tol, kmax=1) {
  #' Find data point nearest to each vertex of network
  #' Assumed 'seg' is sorted in increasing order
  #'         'tp' is increasing within 'seg'
  nX <- length(seg)
  from0 <- from - 1L
  to0   <- to - 1L
  seg0  <- seg - 1L
  #'
  if(kmax == 1) {
    z <- .C("Clinvwhichdist",
            np = as.integer(nX),
            sp = as.integer(seg0),
            tp = as.double(tp),
            nv = as.integer(nv),
            ns = as.integer(ns),
            from = as.integer(from0),
            to   = as.integer(to0),
            seglen = as.double(seglen),
            huge = as.double(huge),
            tol = as.double(tol),
            dist = as.double(numeric(nv)),
            which = as.integer(integer(nv)),
            PACKAGE = "spatstat")
  } else {
    z <- .C("linvknndist",
            kmax = as.integer(kmax), 
            nq = as.integer(nX),
            sq = as.integer(seg0),
            tq = as.double(tp),
            nv = as.integer(nv),
            ns = as.integer(ns),
            from = as.integer(from0),
            to   = as.integer(to0),
            seglen = as.double(seglen),
            huge = as.double(huge),
            tol = as.double(tol),
            dist = as.double(numeric(kmax * nv)),
            which = as.integer(integer(kmax * nv)),
            PACKAGE = "spatstat")
  }
  vnndist <- z$dist
  vnnwhich <- z$which + 1L 
  vnnwhich[vnnwhich == 0] <- NA # possible if network is disconnected
  if(kmax > 1) {
    vnndist <- matrix(vnndist, ncol=kmax, byrow=TRUE)
    vnnwhich <- matrix(vnnwhich, ncol=kmax, byrow=TRUE)
  }
  return(list(vnndist=vnndist, vnnwhich=vnnwhich))
}
