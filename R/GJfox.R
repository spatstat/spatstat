#
#  GJfox.R
#
#  Foxall G-function and J-function
#
#  $Revision: 1.4 $   $Date: 2013/04/25 06:37:43 $
#
Gfox <- function(X, Y, r=NULL, breaks=NULL,
                 correction=c("km", "rs", "han"), ...) {
  stopifnot(is.ppp(X))
  if(!(is.ppp(Y) || is.psp(Y) || is.owin(Y)))
    stop("Y should be an object of class ppp, psp or owin")
  if(!identical(unitname(X), unitname(Y)))
    warning("X and Y are not in the same units")
  # 
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
# ensure compatible windows
  WX <- as.owin(X)
  WY <- as.owin(Y)
  if(!is.subset.owin(WX, WY)) {
    warning("Trimming the window of X to be a subset of the window of Y")
    WX <- intersect.owin(WX, WY)
    X <- X[WX]
  }
# compute distances and censoring distances
  D <- distfun(Y)
  dist <- D(X)
  bdry <- bdist.points(X[WY])
# histogram breakpoints 
  dmax <- max(dist)
  breaks <- handle.r.b.args(r, breaks, WX, NULL, rmaxdefault=dmax)
  rval <- breaks$r
# censoring indicators
  d <- (dist <= bdry)
#  observed distances
  o <- pmin(dist, bdry)
# calculate estimates
  Z <- censtimeCDFest(o, bdry, d, breaks,
                      KM=corx$km,
                      RS=corx$rs,
                      HAN=corx$han,
                      RAW=corx$none,
                      han.denom=if(corx$han) eroded.areas(WX, rval) else NULL)
# relabel
  Z <- rebadge.fv(Z, substitute(Gfox(r), NULL), "Gfox")
  unitname(Z) <- unitname(Y)
  return(Z)
}

Jfox <- function(X, Y, r=NULL, breaks=NULL,
                 correction=c("km", "rs", "han"), ...) {
  H <- Hest(Y, r=r, breaks=breaks, correction=correction, ...)
  G <- Gfox(X, Y, r=H$r, correction=correction, ...)
  # derive J-function
  J <- eval.fv((1-G)/(1-H), dotonly=FALSE)
  # correct calculation of hazard is different
  if("hazard" %in% names(J))
    J$hazard <- G$hazard - H$hazard
  # base labels on 'J' rather than full expression
  attr(J, "labl") <- attr(H, "labl")
  # add column of 1's
  J <- bind.fv(J, data.frame(theo=rep.int(0, nrow(J))), "%s[theo](r)",
               "theoretical value of %s for independence")
  # rename 
  J <- rebadge.fv(J, substitute(Jfox(r), NULL), "Jfox")
  funs <- c("km", "han", "rs", "raw", "theo")
  fvnames(J, ".") <- funs[funs %in% names(J)]
  unitname(J) <- unitname(Y)
  return(J)
}


	

