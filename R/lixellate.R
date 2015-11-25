#'
#'    lixellate.R
#'
#'   Divide each segment of a linear network into several pieces
#' 
#'     $Revision: 1.3 $  $Date: 2015/11/25 04:29:55 $
#'

lixellate <- function(X, ..., nsplit, eps, sparse=TRUE) {
  missn <- missing(nsplit)
  misse <- missing(eps)
  if(missn && misse)
    stop("One of the arguments 'nsplit' or 'eps' must be given")
  if(!missn && !misse)
    stop("The arguments 'nsplit' or 'eps' are incompatible")
  if(!missn)
    stopifnot(is.numeric(nsplit) && all(nsplit >= 0))

  if(is.lpp(X)) {
    rtype <- "lpp"
    np <- npoints(X)
    L <- as.linnet(X)
  } else if(inherits(X, "linnet")) {
    rtype <- "linnet"
    L <- X
    X <- runiflpp(1, L)
    np <- 0
  } else stop("X should be a linnet or lpp object")
  
  if(is.null(sparse))
    sparse <- identical(L$sparse, TRUE)

  from <- L$from
  to <- L$to
  ns <- length(from)

  if(missn) {
    lenfs <- lengths.psp(as.psp(L))
    nsplit <- ceiling(lenfs/eps)
  } else {
    if(length(nsplit) == 1) {
      nsplit <- rep(nsplit, ns)
    } else if(length(nsplit) != ns) {
      stop(paste("nsplit should be a single number,",
                 "or a vector of length equal to the number of segments"))
    }
  }

  sumN <- sum(nsplit)
  sumN1 <- sum(nsplit-1)

  V <- vertices(L)
  nv <- npoints(V)
  xv <- V$x
  yv <- V$y

  coordsX <- coords(X)
  sp <- coordsX$seg
  tp <- coordsX$tp
  ## sort data in increasing order of 'sp'
  oo <- order(sp)
  
  z <- .C("Clixellate",
          ns=as.integer(ns),
          fromcoarse=as.integer(from-1),
          tocoarse = as.integer(to-1),
          fromfine=as.integer(integer(sumN)),
          tofine = as.integer(integer(sumN)),
          nv = as.integer(nv),
          xv = as.double(c(xv, numeric(sumN1))),
          yv = as.double(c(yv, numeric(sumN1))),
          svcoarse = as.integer(integer(nv + sumN1)),
          tvcoarse = as.double(numeric(nv + sumN1)),
          nsplit = as.integer(nsplit),
          np = as.integer(np),
          spcoarse = as.integer(sp[oo]-1),
          tpcoarse = as.double(tp[oo]),
          spfine = as.integer(integer(np)),
          tpfine = as.double(numeric(np)))

  Lfine <- with(z, {
    ii <- seq_len(nv)
    Vnew <- ppp(xv[ii], yv[ii], window=Frame(L), check=FALSE)
    Lfine <- linnet(Vnew, edges=cbind(fromfine,tofine)+1, sparse=sparse)
    marks(Lfine$vertices) <- markcbind(marks(Lfine$vertices),
                                       data.frame(segcoarse=svcoarse+1,
                                                  tpcoarse=tvcoarse))
    Lfine
  })
  if(rtype == "linnet")
    return(Lfine)

  ## put coordinates back in original order
  sp[oo] <- as.integer(z$spfine + 1)
  tp[oo] <- z$tpfine
  coordsX$seg <- sp
  coordsX$tp <- tp
  ## make lpp
  Xfine <- lpp(coordsX, Lfine)
  marks(Xfine) <- marks(X)
  
  return(Xfine)
}

