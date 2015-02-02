#'     quantess.R
#' 
#'     Quantile Tessellation
#'
#'   $Revision: 1.10 $  $Date: 2015/02/01 11:36:57 $

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

nestsplit <- function(X, ...) {
  stopifnot(is.ppp(X))
  flist <- list(...)
  cansplit <- sapply(flist, inherits,
                     what=c("factor", "tess", "owin", "im", "character"))
  splitted <- lapply(flist[cansplit], split, x=X)
  splitters <- lapply(splitted, attr, which="fsplit")
  if(any(!cansplit)) {
    extra <- do.call(MinimalTess, append(list(W=Window(X)), flist[!cansplit]))
    pos <- min(which(!cansplit))
    ns <- length(splitters)
    if(pos > ns) {
      splitters <- append(splitters, list(extra))
    } else {
      before <- splitters[seq_len(pos-1)]
      after  <- splitters[pos:ns]
      splitters <- c(before, list(extra), after)
    }
  }
  ns <- length(splitters)
  if(ns == 0) return(X)
  if(ns == 1) return(split(X, splitters[[1]]))
  if(ns > 2) stop("Nesting depths greater than 2 are not yet implemented")
  names(splitters) <- good.names(names(splitters), paste0("f", 1:ns))
  fax1 <- is.factor(sp1 <- splitters[[1]])
  fax2 <- is.factor(sp2 <- splitters[[2]])
  lev1 <- if(fax1) levels(sp1) else seq_len(sp1$n)
  lev2 <- if(fax2) levels(sp2) else seq_len(sp2$n)
  if(!fax1 && !fax2) {
    ## two tessellations
    marks(sp1) <- factor(lev1, levels=lev1)
    marks(sp2) <- factor(lev2, levels=lev2)
    sp12 <- intersect.tess(sp1, sp2, keepmarks=TRUE)
    pats <- split(X, sp12)
    f1 <- marks(sp12)[,1]
    f2 <- marks(sp12)[,2]
  } else {
    if(fax1 && fax2) {
      ## two grouping factors
      Xsp1 <- split(X, sp1)
      sp2.1 <- split(sp2, sp1)
      ll <- mapply(split, Xsp1, sp2.1, SIMPLIFY=FALSE)
    } else if(fax1 && !fax2) {
      ## grouping factor and tessellation
      Xsp1 <- split(X, sp1)
      ll <- lapply(Xsp1, split, f=sp2)
    } else if(!fax1 && fax2) {
      ## tessellation and grouping factor
      Xsp1 <- split(X, sp1)
      sp2.1 <- split(sp2, attr(Xsp1, "fgroup"))
      ll <- mapply(split, Xsp1, sp2.1, SIMPLIFY=FALSE)
    }
    neach <- sapply(ll, length)
    f1 <- rep(factor(lev1, levels=lev1), neach)
    f2 <- rep(factor(lev2, levels=lev2), length(Xsp1))
    pats <- do.call(c, unname(ll))
  }
  h <- hyperframe(pts=pats, f1=f1, f2=f2)
  names(h)[2:3] <- names(splitters)
  return(h)
}
