#'  tapplysum.R
#'
#'  Faster replacement for tapply(..., FUN=sum)
#'
#'  Adrian Baddeley and Tilman Davies
#' 
#'  $Revision: 1.11 $  $Date: 2016/12/12 09:07:06 $

tapplysum <- function(x, flist, do.names=FALSE, na.rm=TRUE) {
  stopifnot(is.numeric(x))
  stopifnot(is.list(flist))
  stopifnot(all(lengths(flist) == length(x)))
  stopifnot(all(sapply(flist, is.factor)))
  nfac <- length(flist)
  goodx <- is.finite(x)
  if(na.rm) goodx <- goodx | is.na(x)
  if(!(nfac %in% 2:3) || !all(goodx)) {
    y <- tapply(x, flist, sum)
    y[is.na(y)] <- 0
    return(y)
  }
  ifac <- flist[[1L]]
  jfac <- flist[[2L]]
  mi <- length(levels(ifac))
  mj <- length(levels(jfac))
  ii <- as.integer(ifac)
  jj <- as.integer(jfac)
  if(nfac == 3) {
    kfac <- flist[[3L]]
    mk <- length(levels(kfac))
    kk <- as.integer(kfac)
  }
  #' remove NA's
  if(nfac == 2) {
    if(anyNA(x) || anyNA(ii) || anyNA(jj)) {
      ok <- !(is.na(x) | is.na(ii) | is.na(jj))
      ii <- ii[ok]
      jj <- jj[ok]
      x <- x[ok]
    }
  } else {
    if(anyNA(x) || anyNA(ii) || anyNA(jj) || anyNA(kk)) {
      ok <- !(is.na(x) | is.na(ii) | is.na(jj) | is.na(kk))
      ii <- ii[ok]
      jj <- jj[ok]
      kk <- kk[ok]
      x <- x[ok]
    }
  }
  n <- length(ii)
  #' 
  if(nfac == 2) {
    result <- matrix(0, mi, mj)
    if(n > 0) {
      oo <- order(ii, jj)
      zz <- .C("ply2sum",
               nin = as.integer(n),
               xin = as.double(x[oo]),
               iin = as.integer(ii[oo]),
               jin = as.integer(jj[oo]),
               nout = as.integer(integer(1L)),
               xout = as.double(numeric(n)),
               iout = as.integer(integer(n)),
               jout = as.integer(integer(n)))
      nout <- zz$nout
      if(nout > 0) {
        ijout <- cbind(zz$iout, zz$jout)[1:nout,,drop=FALSE]
        xout  <- zz$xout[1:nout]
        result[ijout] <- xout
      }
    }
  } else {
    result <- array(0, dim=c(mi, mj, mk))
    if(n > 0) {
      oo <- order(ii, jj, kk)
      zz <- .C("ply3sum",
               nin = as.integer(n),
               xin = as.double(x[oo]),
               iin = as.integer(ii[oo]),
               jin = as.integer(jj[oo]),
               kin = as.integer(kk[oo]),
               nout = as.integer(integer(1L)),
               xout = as.double(numeric(n)),
               iout = as.integer(integer(n)),
               jout = as.integer(integer(n)),
               kout = as.integer(integer(n)))
      nout <- zz$nout
      if(nout > 0) {
        ijkout <- cbind(zz$iout, zz$jout, zz$kout)[1:nout,,drop=FALSE]
        xout  <- zz$xout[1:nout]
        result[ijkout] <- xout
      }
    }
  }
  if(do.names) 
    dimnames(result) <- lapply(flist, levels)
  return(result)
}


                       
           
