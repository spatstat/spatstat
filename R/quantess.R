#'     quantess.R
#' 
#'     Quantile Tessellation
#'
#'   $Revision: 1.18 $  $Date: 2019/03/21 02:56:53 $

quantess <- function(M, Z, n, ...) {
  UseMethod("quantess")
}

quantess.owin <- function(M, Z, n, ..., type=2, origin=c(0,0), eps=NULL) {
  W <- as.owin(M)
  B <- boundingbox(W)
  tcross <- MinimalTess(W, ...)
  force(n)
  if(!is.character(Z)) {
    Zim <- as.im(Z, W, eps=eps)
    Zrange <- range(Zim)
  } else {
    Z <- match.arg(Z, c("x", "y", "rad", "ang"))
    if(Z %in% c("x", "y") && is.rectangle(W)) {
      out <- switch(Z,
                    x={ quadrats(W, nx=n, ny=1) },
                    y={ quadrats(W, nx=1, ny=n) })
      if(!is.null(tcross)) out <- intersect.tess(out, tcross)
      return(out)
    }
    a <- qtPrepareCoordinate(Z, W, origin)
    Zfun <- a$Zfun
    Zrange <- a$Zrange
    Zim <- as.im(Zfun, W, eps=eps)
  }
  qZ <- quantile(Zim, probs=(0:n)/n, type=type)
  qZ[1] <- min(qZ[1], Zrange[1])
  qZ[n+1] <- max(qZ[n+1], Zrange[2])
  if(is.polygonal(W) && is.character(Z)) {
    R <- Frame(W)
    strips <- switch(Z,
                     x = tess(xgrid=qZ, ygrid=R$yrange),
                     y = tess(xgrid=R$xrange, ygrid=qZ),
                     rad = polartess(B, radii=qZ, origin=origin),
                     ang = polartess(B, angles=qZ, origin=origin))
    out <- intersect.tess(strips, tess(tiles=list(W)))
    qzz <- signif(qZ, 3)
    tilenames(out) <- paste0("[", qzz[1:n], ",",
                             qzz[-1], c(rep(")", n-1), "]"))
  } else {
    ZC <- cut(Zim, breaks=qZ, include.lowest=TRUE, right=FALSE)
    out <- tess(image=ZC)
  }
  if(!is.null(tcross)) out <- intersect.tess(out, tcross)
  return(out)
}

qtPrepareCoordinate <- function(covname, W, origin=c(0,0)) {
  switch(covname,
         x={
           Zfun <- function(x,y){x}
           Zrange <- boundingbox(W)$xrange
         },
         y={
           Zfun <- function(x,y){y}
           Zrange <- boundingbox(W)$yrange
         },
         rad={
           origin <- interpretAsOrigin(origin, W)
           Zfun <- function(x,y) { sqrt((x-origin[1])^2+(y-origin[2])^2) }
           V <- vertices(W)
           Zrange <- range(Zfun(V$x, V$y))
         },
         ang={
           origin <- interpretAsOrigin(origin, W)
           Zstart <- 0
           Zfun <- function(x,y) {
             angle <- atan2(y-origin[2], x-origin[1]) %% (2*pi)
             if(Zstart < 0) {
               negangle <- angle - 2*pi
               angle <- ifelse(negangle >= Zstart, negangle, angle)
             }
             return(angle)
           }
           S <- as.data.frame(edges(W))
           a <- Zfun(S[,"x0"], S[,"y0"])
           b <- Zfun(S[,"x1"], S[,"y1"])
           bmina <- b - a
           swap <- (bmina > pi) | (bmina < 0 & bmina > -pi)
           arcs <- cbind(ifelse(swap, b, a), ifelse(swap, a, b))
           arcs <- lapply(apply(arcs, 1, list), unlist)
           Zunion <- circunion(arcs)
           Zrange <- c(Zunion[[1]][1], Zunion[[length(Zunion)]][2])
           if(diff(Zrange) < 0) {
             #' first interval straddles the positive x-axis
             Zstart <- Zrange[1] <- Zrange[1] - 2*pi
           } 
         })
  return(list(Zrange=Zrange, Zfun=Zfun))
}

quantess.ppp <- function(M, Z, n, ..., type=2, origin=c(0,0), eps=NULL) {
  W <- as.owin(M)
  B <- boundingbox(W)
  tcross <- MinimalTess(W, ...)
  force(n)
  if(!is.character(Z)) {
    Zim <- as.im(Z, W, eps=eps)
    ZM <- if(is.function(Z)) Z(M$x, M$y) else Zim[M]
    Zrange <- range(range(Zim), ZM)
  } else {
    Z <- match.arg(Z, c("x", "y", "rad", "ang"))
    if(Z %in% c("x", "y") && is.rectangle(W)) {
      switch(Z,
             x={
               qx <- quantile(M$x, probs=(1:(n-1))/n, type=type)
               qx <- c(W$xrange[1], qx, W$xrange[2])
               out <- tess(xgrid=qx, ygrid=W$yrange)
             },
             y={
               qy <- quantile(M$y, probs=(1:(n-1))/n, type=type)
               qy <- c(W$yrange[1], qy, W$yrange[2])
               out <- tess(xgrid=W$xrange, ygrid=qy)
             })
      if(!is.null(tcross)) out <- intersect.tess(out, tcross)
      return(out)
    }
    a <- qtPrepareCoordinate(Z, W, origin)
    Zrange <- a$Zrange
    Zfun <- a$Zfun
    ZM <- Zfun(M$x, M$y)
    Zrange <- range(Zrange, range(ZM))
    Zim <- as.im(Zfun, W, eps=eps)
  } 
  qZ <- quantile(Zim, probs=(0:n)/n, type=type)
  qZ[1] <- min(qZ[1], Zrange[1])
  qZ[n+1] <- max(qZ[n+1], Zrange[2])
  if(is.polygonal(W) && is.character(Z)) {
    R <- Frame(W)
    strips <- switch(Z,
                     x = tess(xgrid=qZ, ygrid=R$yrange),
                     y = tess(xgrid=R$xrange, ygrid=qZ),
                     rad = polartess(B, radii=qZ, origin=origin),
                     ang = polartess(B, angles=qZ, origin=origin))
    out <- intersect.tess(strips, tess(tiles=list(W)))
    qzz <- signif(qZ, 3)
    tilenames(out) <- paste0("[", qzz[1:n], ",",
                             qzz[-1], c(rep(")", n-1), "]"))
  } else {
    ZC <- cut(Zim, breaks=qZ, include.lowest=TRUE)
    out <- tess(image=ZC)
  }
  if(!is.null(tcross)) out <- intersect.tess(out, tcross)
  return(out)
}

quantess.im <- function(M, Z, n, ..., type=2, origin=c(0,0)) {
  W <- Window(M)
  tcross <- MinimalTess(W, ...)
  force(n)
  if(!(type %in% c(1,2)))
    stop("Only quantiles of type 1 and 2 are implemented for quantess.im")
  if(is.character(Z)) {
    Z <- match.arg(Z, c("x", "y", "rad", "ang"))
    a <- qtPrepareCoordinate(Z, W, origin)
    Z <- a$Zfun
    Zrange <- a$Zrange
  } else Zrange <- NULL
  MZ <- harmonise(M=M, Z=Z)
  M <- MZ$M[W, drop=FALSE]
  Z <- MZ$Z[W, drop=FALSE]
  Zrange <- range(c(range(Z), Zrange))
  Fun <- ewcdf(Z[], weights=M[]/sum(M[]))
  qZ <- quantile(Fun, probs=(1:(n-1))/n, type=type)
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
  ## find the minimal tessellation of W consistent with the arguments
  argh <- list(...)
  v <- NULL
  if(length(argh)) {
    nama <- names(argh)
    known <- union(names(formals(quadrats)),
                   names(formals(tess)))
    recognised <- !is.na(match(nama, known))
    if(any(recognised)) {
      if(any(c("nx", "ny") %in% nama)) {
        v <- do.call(quadrats,
                     resolve.defaults(list(X=W),
                                      argh[recognised],
                                      list(nx=1, ny=1)))
      } else if(any(c("xbreaks", "ybreaks") %in% nama)) {
        v <- do.call(quadrats,
                     resolve.defaults(list(X=W),
                                      argh[recognised],
                                      list(xbreaks=W$xrange,
                                           ybreaks=W$yrange)))
      } else {
        v <- do.call(tess,
                     resolve.defaults(argh[recognised],
                                      list(window=W,
                                           keepempty=TRUE)))
      }
    }
  } 
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
    neach <- lengths(ll)
    f1 <- rep(factor(lev1, levels=lev1), neach)
    f2 <- rep(factor(lev2, levels=lev2), length(Xsp1))
    pats <- do.call(c, unname(ll))
  }
  h <- hyperframe(pts=pats, f1=f1, f2=f2)
  names(h)[2:3] <- names(splitters)
  return(h)
}
