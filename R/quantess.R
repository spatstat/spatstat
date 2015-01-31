#'     quantess.R
#' 
#'     Quantile Tessellation
#'
#'   $Revision: 1.7 $  $Date: 2015/01/31 14:21:13 $

quantess <- function(M, Z, n, ...) {
  UseMethod("quantess")
}

quantess.owin <- function(M, Z, n, ...) {
  W <- as.owin(M)
  tcross <- MinimalTess(W, ...)
  if(!is.character(Z)) {
    Zim <- as.im(Z, W)
    Zrange <- range(Zim)
  } else {
    if(!(Z %in% c("x", "y")))
      stop(paste("Unrecognised covariate", dQuote(Z)))
    if(is.rectangle(W)) {
      out <- switch(Z,
                    x={ quadrats(W, nx=n, ny=1) },
                    y={ quadrats(W, nx=1, ny=n) })
      if(!is.null(tcross)) out <- intersect.tess(out, tcross)
      return(out)
    }
    switch(Z,
           x={
             Zfun <- function(x,y){x}
             Zrange <- boundingbox(W)$xrange
           },
           y={
             Zfun <- function(x,y){y}
             Zrange <- boundingbox(W)$yrange
           })
    Zim <- as.im(Zfun, W)
  }
  qZ <- quantile(Zim, probs=(1:(n-1))/n)
  qZ <- c(Zrange[1], qZ, Zrange[2])
  if(is.polygonal(W) && is.character(Z)) {
    R <- Frame(W)
    strips <- switch(Z,
                     x = tess(xgrid=qZ, ygrid=R$yrange),
                     y = tess(xgrid=R$xrange, ygrid=qZ))
    out <- intersect.tess(strips, tess(tiles=list(W)))
  } else {
    ZC <- cut(Z, breaks=qZ, include.lowest=TRUE)
    out <- tess(image=ZC)
  }
  qzz <- signif(qZ, 3)
  tilenames(out) <- paste0("[", qzz[1:(n-1)], ",",
                           qzz[-1], c(rep(")", n-1), "]"))
  if(!is.null(tcross)) out <- intersect.tess(out, tcross)
  return(out)
}

quantess.ppp <- function(M, Z, n, ...) {
  W <- as.owin(M)
  tcross <- MinimalTess(W, ...)
  if(!is.character(Z)) {
    Zim <- as.im(Z, W)
    ZM <- if(is.function(Z)) Z(M$x, M$y) else Zim[M]
    Zrange <- range(range(Zim), ZM)
  } else {
    if(!(Z %in% c("x", "y")))
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
    switch(Z,
           x={
             Zfun <- function(x,y){x}
             ZM <- M$x
             Zrange <- boundingbox(W)$xrange
           },
           y={
             Zfun <- function(x,y){y}
             ZM <- M$y
             Zrange <- boundingbox(W)$yrange
           })
    Zim <- as.im(Zfun, W)
  } 
  qZ <- quantile(ZM, probs=(1:(n-1))/n)
  qZ <- c(Zrange[1], qZ, Zrange[2])
  if(is.polygonal(W) && is.character(Z)) {
    R <- Frame(W)
    strips <- switch(Z,
                     x = tess(xgrid=qZ, ygrid=R$yrange),
                     y = tess(xgrid=R$xrange, ygrid=qZ))
    out <- intersect.tess(strips, tess(tiles=list(W)))
  } else {
    ZC <- cut(Zim, breaks=qZ, include.lowest=TRUE)
    out <- tess(image=ZC)
  }
  qzz <- signif(qZ, 3)
  tilenames(out) <- paste0("[", qzz[1:(n-1)], ",",
                           qzz[-1], c(rep(")", n-1), "]"))
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
  Zrange <- range(Z)
  Fun <- ewcdf(Z[], weights=M[]/sum(M[]))
  qZ <- quantile(Fun, probs=(1:(n-1))/n)
  qZ <- c(Zrange[1], qZ, Zrange[2])
  ZC <- cut(Z, breaks=qZ, include.lowest=TRUE)
  out <- tess(image=ZC)
  qzz <- signif(qZ, 3)
  tilenames(out) <- paste0("[", qzz[1:(n-1)], ",",
                           qzz[-1], c(rep(")", n-1), "]"))
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

quantsplit <- function(M, Z, n, ...) {
  f <- quantess(M, Z, n, ...)
  out <- if(is.owin(M)) f else split(M, f)
  return(out)
}
