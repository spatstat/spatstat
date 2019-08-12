#
#  GJfox.R
#
#  Foxall G-function and J-function
#
#  $Revision: 1.10 $   $Date: 2019/08/12 10:10:00 $
#
Gfox <- function(X, Y, r=NULL, breaks=NULL,
                 correction=c("km", "rs", "han"),
                 W=NULL, ...) {
  stopifnot(is.ppp(X))
  #' validate and resolve windows
  a <- resolve.foxall.window(X, Y, W)
  X <- a$X
  Y <- a$Y
  W <- a$W
  #' 
  if(is.null(correction))
    correction <- c("rs", "km", "cs")
  correction <- pickoption("correction", correction,
                           c(none="none",
                             raw="none",
                             border="rs",
                             rs="rs",
                             KM="km",
                             km="km",
                             Kaplan="km",
                             han="han",
                             Hanisch="han",
                             best="km"),
                           multi=TRUE)
  corxtable <- c("km", "rs", "han", "none") 
  corx <- as.list(corxtable %in% correction)
  names(corx) <- corxtable
  ## compute distances and censoring distances
  D <- distfun(Y)
  dist <- D(X)
  bdry <- bdist.points(X[W]) # sic
  ## histogram breakpoints 
  dmax <- max(dist)
  breaks <- handle.r.b.args(r, breaks, Window(X), NULL, rmaxdefault=dmax)
  rval <- breaks$r
  ## censoring indicators
  d <- (dist <= bdry)
  ## observed distances
  o <- pmin.int(dist, bdry)
  ## calculate estimates
  Z <- censtimeCDFest(o, bdry, d, breaks,
                      KM=corx$km,
                      RS=corx$rs,
                      HAN=corx$han,
                      RAW=corx$none,
                      han.denom=if(corx$han) eroded.areas(Window(X), rval) else NULL,
                      tt=dist)
  ## relabel
  Z <- rebadge.fv(Z, quote(G[fox](r)), c("G", "fox"))
  unitname(Z) <- unitname(Y)
  return(Z)
}

Jfox <- function(X, Y, r=NULL, breaks=NULL,
                 correction=c("km", "rs", "han"), W=NULL, ...) {
  ## validate and resolve windows
  a <- resolve.foxall.window(X, Y, W)
  X <- a$X
  Y <- a$Y
  W <- a$W
  ##  process
  H <- Hest(Y, r=r, breaks=breaks, correction=correction, ..., W=W)
  G <- Gfox(X, Y, r=H$r, correction=correction, ..., W=W)
  ## derive J-function
  J <- eval.fv((1-G)/(1-H), dotonly=FALSE)
  ## correct calculation of hazard is different
  if("hazard" %in% names(J))
    J$hazard <- G$hazard - H$hazard
  ## base labels on 'J' rather than full expression
  attr(J, "labl") <- attr(H, "labl")
  ## add column of 1's
  J <- bind.fv(J, data.frame(theo=rep.int(1, nrow(J))), "%s[theo](r)",
               "theoretical value of %s for independence")
  ## rename 
  J <- rebadge.fv(J, quote(J[fox](r)), c("J", "fox"))
  funs <- c("km", "han", "rs", "raw", "theo")
  fvnames(J, ".") <- funs[funs %in% names(J)]
  unitname(J) <- unitname(Y)
  return(J)
}

resolve.foxall.window <- function(X, Y, W=NULL) {
  if(!(is.ppp(Y) || is.psp(Y) || is.owin(Y) || is.im(Y)))
    stop("Y should be an object of class ppp, psp, owin or im")
  if(is.im(Y) && !is.logical(ZeroValue(Y)))
    stop("When Y is an image, its pixel values should be logical values")
  if(!identical(unitname(X), unitname(Y)))
    warning("X and Y are not in the same units")
  ## default window based on Y
  if(is.ppp(Y) || is.psp(Y)) {
    W0 <- Window(Y)
    W0describe <- "the observation window of Y"
  } else if(is.owin(Y)) {
    W0 <- Frame(Y)
    W0describe <- "the Frame of Y"
  } else if(is.im(Y)) {
    W0 <- Window(Y)
    W0describe <- "the observation window of Y"
    Y <- solutionset(Y)
  } else stop("Y should be an object of class ppp, psp, owin or im")
  ## actual window used for estimation
  if(!is.null(W)) {
    stopifnot(is.owin(W))
    if(!is.subset.owin(W, W0))
      stop(paste("W is not a subset of", W0describe))
    Wdescribe <- "W"
  } else {
    W <- W0
    Wdescribe <- W0describe
  }
  ## ensure compatible windows
  WX <- Window(X)
  if(!is.subset.owin(WX, W)) {
    warning(paste("Trimming the window of X to be a subset of", Wdescribe))
    WX <- intersect.owin(WX, W)
    if(area.owin(WX) == 0) stop("Trimmed window has zero area")
    X <- X[WX]
    if(npoints(X) == 0) stop("No points remaining after trimming window")
  }
  return(list(X=X, Y=Y, W=W))
}
