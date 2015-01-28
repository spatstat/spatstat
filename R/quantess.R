#'     quantess.R
#' 
#'     Quantile Tessellation
#'
#'   $Revision: 1.2 $  $Date: 2015/01/28 12:04:44 $

quantess <- function(M, Z, n, ...) {
  UseMethod("quantess")
}

quantess.owin <- function(M, Z, n, ...) {
  W <- as.owin(M)
  tcross <- MinimalTess(W, ...)
  if(is.character(Z)) {
    if(!any(c("x", "y") == Z))
      stop(paste("Unrecognised covariate", dQuote(Z)))
    if(is.rectangle(W)) {
      out <- switch(Z,
                    x={ quadrats(W, nx=n, ny=1) },
                    y={ quadrats(W, nx=1, ny=n) })
      if(!is.null(tcross)) out <- intersect.tess(out, tcross)
      return(out)
    }
    Z <- switch(Z,
                x=function(x,y){x},
                y=function(x,y){y})
  }
  Z <- as.im(Z, W)
  qZ <- quantile(Z, probs=(1:(n-1))/n)
  qZ <- c(min(Z), qZ, max(Z))
  ZC <- cut(Z, breaks=qZ, include.lowest=TRUE)
  out <- tess(image=ZC)
  if(!is.null(tcross)) out <- intersect.tess(out, tcross)
  return(out)
}

quantess.ppp <- function(M, Z, n, ...) {
  W <- as.owin(M)
  tcross <- MinimalTess(W, ...)
  if(is.character(Z)) {
    if(!any(c("x", "y") == Z))
      stop(paste("Unrecognised covariate", dQuote(Z)))
    if(is.rectangle(W)) {
      switch(Z,
             x={
               qx <- quantile(M$x, probs=(1:(n-1))/n)
               qx <- c(W$xrange[1], qx, W$xrange[2])
               out <- tess(xgrid=qx, ygrid=W$yrange)
             },
             y={
               qy <- quantile(M$y, probs=(1:(n-1))/n)
               qy <- c(W$yrange[1], qy, W$yrange[2])
               out <- tess(xgrid=W$xrange, ygrid=qy)
             })
      if(!is.null(tcross)) out <- intersect.tess(out, tcross)
      return(out)
    }
    Z <- switch(Z,
                x=function(x,y){x},
                y=function(x,y){y})
  }
  Z <- as.im(Z, W)
  ZM <- Z[M]
  qZ <- quantile(ZM, probs=(1:(n-1))/n)
  qZ <- c(min(ZM), qZ, max(ZM))
  ZC <- cut(Z, breaks=qZ, include.lowest=TRUE)
  out <- tess(image=ZC)
  if(!is.null(tcross)) out <- intersect.tess(out, tcross)
  return(out)
}

quantess.im <- function(M, Z, n, ...) {
  W <- Window(M)
  tcross <- MinimalTess(W, ...)
  if(is.character(Z)) 
    Z <- switch(Z,
                x=function(x,y){x},
                y=function(x,y){y},
                stop(paste("Unrecognised covariate", dQuote(Z))))
  MZ <- harmonise(M=M, Z=Z)
  M <- MZ$M[W, drop=FALSE]
  Z <- MZ$Z[W, drop=FALSE]
  Fun <- ewcdf(Z[], weights=M[]/sum(M[]))
  qZ <- quantile(Fun, probs=(1:(n-1))/n)
  qZ <- c(min(Z), qZ, max(Z))
  ZC <- cut(Z, breaks=qZ, include.lowest=TRUE)
  out <- tess(image=ZC)
  if(!is.null(tcross)) out <- intersect.tess(out, tcross)
  return(out)
}

MinimalTess <- function(W, ...) {
  # find the minimal tessellation of W consistent with the arguments 
  argh <- list(...)
  if(length(argh) == 0) return(NULL)
  nama <- names(argh)
  if(any(c("nx", "ny") %in% nama)) {
    fun <- quadrats
    dflt <- list(nx=1, ny=1)
  } else if(any(c("xbreaks", "ybreaks") %in% nama)) {
    fun <- quadrats
    dflt <- list(xbreaks=W$xrange, ybreaks=W$yrange)
  } else {
    fun <- tess
    dflt <- list(window=W, keepempty=TRUE)
  }
  v <- do.call(fun, resolve.defaults(list(W), argh, dflt))
  return(v)
}
