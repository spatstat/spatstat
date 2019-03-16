#'
#'   polartess.R
#'
#'   Tessellation using polar coordinates
#' 
#'   $Revision: 1.4 $  $Date: 2019/03/16 05:36:40 $

polartess <- function(W, ..., nradial=NULL, nangular=NULL,
                      radii=NULL, angles=NULL, origin=NULL,
                      sep="x") {
  trap.extra.arguments(...)
  W <- as.owin(W)
  if(!is.null(origin)) {
    origin <- interpretAsOrigin(origin, W)
    W <- shift(W, -origin)
  }

  V <- vertices(Frame(W))
  rmax <- sqrt(max(V$x^2 + V$y^2))

  if(!is.null(radii)) {
    if(!is.null(nradial))
      warning("nradial ignored because radii were specified")
    radii <- as.numeric(radii)
    stopifnot(length(radii) >= 2)
    stopifnot(all(radii >= 0))
    if(sum(is.infinite(radii)) > 1 || !all(diff(radii) > 0))
      stop("radii should be increasing")
    radnames <- paste(signif(radii, 4))
    radii[is.infinite(radii)] <- 1.01 * rmax
    rmax <- max(radii)
    nradial <- length(radii) - 1L
  } else if(!is.null(nradial)) {
    check.1.integer(nradial)
    radii <- seq(0, rmax, length.out=nradial+1L)
    radnames <- paste(signif(radii, 4))
  }
  nradii <- length(radii)
  
  if(!is.null(angles)) {
    if(!is.null(nangular))
      warning("nangular ignored because angles were specified")
    angles <- as.numeric(angles)
    stopifnot(length(angles) >= 2)
    if(!all(diff(angles) > 0))
      stop("angles should be increasing")
    if(diff(range(angles)) > 2 * pi + .Machine$double.eps)
      stop("The range of angles must not exceed 2 * pi")
    nangular <- length(angles) - 1L
  } else if(!is.null(nangular)) {
    check.1.integer(nangular)
    angles <- seq(0, 2*pi, length.out=nangular+1L)
  }
  nangles <- length(angles)
  
  #' build tessellations
  result <- as.tess(W)
  DD <- Dmax <- disc(rmax)
  
  if(!is.null(radii)) {
    rmin <- radii[1]
    if(rmin > 0) DD <- setminus.owin(DD, disc(rmin))
    Dlist <- lapply(radii[radii > 0], disc)
    if(rmin == 0) Dlist <- append(list(NULL), Dlist)
    Tlist <- list()
    for(i in 1:nradial) 
      Tlist <- append(Tlist, list(setminus.owin(Dlist[[i+1]], Dlist[[i]])))
    names(Tlist) <- paste0("[", radnames[-nradii], ", ", radnames[-1L],
                           c(rep(")", nradial-1L), "]"))
    Rtess <- tess(tiles=Tlist, window=DD)
    result <- intersect.tess(result, Rtess, sep=sep)
  } 
  
  if(!is.null(angles)) {
    Tlist <- list()
    aa <- seq(min(angles), max(angles), length.out=256)
    aa <- sort(c(aa, angles))
    xx <- rmax * cos(aa)
    yy <- rmax * sin(aa)
    for(i in 1:nangular) {
      jj <- (aa >= angles[i]) & (aa <= angles[i+1L])
      Tlist[[i]] <- owin(poly=list(x=c(0, xx[jj]),
                                   y=c(0, yy[jj])))
    }
    angnames <- lapply(angles/pi, simplenumber, unit="pi", multiply="")
    unknown <- sapply(angnames, is.null)
    angnames[unknown] <- paste(signif((angles/pi)[unknown], 4), "pi")
    angnames <- unlist(angnames)
    names(Tlist) <- paste0("[", angnames[-nangles], ", ", angnames[-1L],
                           c(rep(")", nangular-1L), "]"))
    gap <- abs(1 - diff(range(angles))/(2*pi))
    DDD <- if(gap < 0.01) Dmax else owin(poly=list(x=c(0, xx), y=c(0,yy)))
    Atess <- tess(tiles=Tlist, window=DDD)
    result <- intersect.tess(result, Atess, sep=sep)
  }
  if(!is.null(origin))
    result <- shift(result, vec=origin)
  return(result)
}
