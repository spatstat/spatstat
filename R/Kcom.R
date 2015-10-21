#
#  Kcom.R
#
#   model compensated K-function
#
# $Revision: 1.14 $ $Date: 2015/10/21 09:06:57 $
#

Kcom <- local({

  Kcom <- function(object, r=NULL, breaks=NULL, ..., 
                   correction=c("border", "isotropic", "translate"),
                   conditional=!is.poisson(object),
                   restrict=FALSE,
                   model=NULL, 
                   trend=~1, interaction=Poisson(), rbord=reach(interaction),
                   compute.var=TRUE,
                   truecoef=NULL, hi.res=NULL) {
  if(inherits(object, "ppm")) {
    fit <- object
  } else if(is.ppp(object) || inherits(object, "quad")) {
    if(is.ppp(object)) object <- quadscheme(object, ...)
    if(!is.null(model)) {
      fit <- update(model, Q=object, forcefit=TRUE)
    } else {
      fit <- ppm(object, trend=trend, interaction=interaction, rbord=rbord,
                 forcefit=TRUE)
    }
  } else 
    stop("object should be a fitted point process model or a point pattern")
  
  if(missing(conditional) || is.null(conditional))
    conditional <- !is.poisson(fit)

#  rfixed <- !is.null(r) || !is.null(breaks)
  
  # Extract data and window
  Q <- quad.ppm(fit, drop=FALSE)
  X <- data.ppm(fit)
  Win <- X$window

  # selection of edge corrections
  correction.given <- !missing(correction) && !is.null(correction)
  correction <- pickoption("correction", correction,
                           c(none="none",
                             border="border",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             ripley="isotropic",
                             trans="translation",
                             translate="translation",
                             translation="translation",
                             best="best"),
                           multi=TRUE)
  correction <- implemented.for.K(correction, Win$type, correction.given)

  opt <- list(bord = any(correction == "border"),
              tran = any(correction == "translation"),
              ripl = any(correction == "isotropic"))
  if(sum(unlist(opt)) == 0)
    stop("No corrections selected")
  
  # edge correction algorithm
  algo <- if(!conditional) "classical" else
          if(restrict) "restricted" else "reweighted"

  # conditioning on border region?
  if(!conditional) {
    Wfree <- Win
  } else {
    rbord <- fit$rbord
    Wfree <- erosion(Win, rbord)
    if(restrict) {
      retain <- inside.owin(union.quad(Q), , Wfree)
      # Throw away boundary data
      Q <- Q[Wfree]
      X <- X[Wfree]
      Win <- Wfree
    }
  }
  
  # Extract quadrature info
  U <- union.quad(Q)
  Z <- is.data(Q) # indicator data/dummy
  E <- equalsfun.quad(Q)
  WQ <- w.quad(Q)  # quadrature weights

  # quadrature points used 
  USED <- if(algo == "reweighted") (bdist.points(U) > rbord) else rep.int(TRUE, U$n)

  # basic statistics
  npts <- npoints(X)
  areaW <- area(Win)
  lambda <- npts/areaW
  lambda2 <- npts * (npts - 1)/(areaW^2)

  # adjustments to account for restricted domain of pseudolikelihood
  if(algo == "reweighted") {
    npts.used <- sum(Z & USED)
    area.used <- sum(WQ[USED])
#    lambda.used <- npts.used/area.used
#    lambda2.used <- npts.used * (npts.used - 1)/(area.used^2)
  } else {
    npts.used <- npts
    area.used <- areaW
#    lambda.used <- lambda
#    lambda2.used <- lambda2
  }
  
  # 'r' values
  rmaxdefault <- rmax.rule("K", if(restrict) Wfree else Win, npts/areaW)
  breaks <- handle.r.b.args(r, breaks, Wfree, rmaxdefault=rmaxdefault)
  r <- breaks$r
#  nr <- length(r)
  rmax <- breaks$max

  
  # recommended range of r values
  alim <- c(0, min(rmax, rmaxdefault))
        
  # this will be the output data frame
  K <- data.frame(r=r, pois=pi * r^2)
  desc <- c("distance argument r", "expected %s for CSR")
  K <- fv(K, "r", substitute(K(r), NULL),
            "pois", , alim, c("r","%s[pois](r)"), desc, fname="K")

  ############### start computing ##################

  # residuals
  resid <- residuals(fit, type="raw",drop=FALSE,
                    new.coef=truecoef, quad=hi.res)
  resval  <- with(resid, "increment")
  rescts  <- with(resid, "continuous")
  if(restrict) {
    # keep only data inside Wfree
    resval <- resval[retain]
    rescts <- rescts[retain]
  }
  
  # close pairs of points
  # (quadrature point to data point)
  clos <- crosspairs(U, X, rmax, what="ijd")
  dIJ <- clos$d
  I   <- clos$i
  J   <- clos$j
  UI <- U[I]
  XJ <- X[J]
  EIJ <- E(I, J) # TRUE if points are identical, U[I[k]] == X[J[k]] 
  ZI <- Z[I]     # TRUE if U[I[k]] is a data point
  DD <- ZI & !EIJ  # TRUE for pairs of distinct data points only
#  nDD <- sum(DD)

  # determine whether a quadrature point will be used in integral
  okI <- USED[I]
  
  if(spatstat.options("Kcom.remove.zeroes"))
    okI <- okI & !EIJ
  
  # residual weights
#  wIJ <- ifelseXY(EIJ, rescts[I], resval[I])
  # absolute weight for continuous integrals
  wc   <- -rescts
  wcIJ <- -rescts[I]

  ####################################################
  
  if(opt$bord) {
    # border method
    # Compute distances to boundary
    # (in restricted case, the window of U has been adjusted)
    b <- bdist.points(U)
    bI <- b[I]
    # reduced sample for K(r) of data only
    RSX <- Kount(dIJ[DD & okI], bI[DD & okI], b[Z & USED], breaks)
#    Kb <- RSX$numerator/(lambda.used * RSX$denom.count)
    Kb <- RSX$numerator/(lambda * RSX$denom.count)
    K <- bind.fv(K, data.frame(border=Kb), "hat(%s)[bord](r)",
                 nzpaste(algo,
                         "border-corrected nonparametric estimate of %s"),
                 "border")
    # reduced sample for adjustment integral
    RSD <- Kwtsum(dIJ[okI], bI[okI], wcIJ[okI],
                  b[Z & USED], rep.int(1, npts.used), breaks)
#    lambdaU <- (npts.used + 1)/area.used
    lambdaU <- (npts + 1)/areaW
    Kb <- RSD$numerator/((RSD$denominator + 1) * lambdaU)

    K <- bind.fv(K, data.frame(bcom=Kb), "bold(C)~hat(%s)[bord](r)",
                 nzpaste("model compensator of",
                         algo, "border-corrected %s"),
                 "border")
  }
  if(opt$tran) {
    # translation correction
    edgewt <- switch(algo,
                     classical  = edge.Trans(UI, XJ, paired=TRUE),
                     restricted = edge.Trans(UI, XJ, paired=TRUE),
                     reweighted = edge.Trans.modif(UI, XJ, Win, Wfree,
                       paired=TRUE))
    wh   <- whist(dIJ[okI], breaks$val, (edgewt * wcIJ)[okI])
    whDD <- whist(dIJ[DD & okI], breaks$val, edgewt[DD & okI])    
    Ktrans <- cumsum(whDD)/(lambda2 * area.used)
    Ktrans[r >= rmax] <- NA
    K <- bind.fv(K, data.frame(trans=Ktrans), "hat(%s)[trans](r)",
                 nzpaste(algo,
                         "translation-corrected nonparametric estimate of %s"),
                 "trans")
#    lambda2U <- (npts.used + 1) * npts.used/(area.used^2)
    lambda2U <- (npts + 1) * npts/(areaW^2)
    Ktrans <- cumsum(wh)/(lambda2U * area.used)
    Ktrans[r >= rmax] <- NA
    K <- bind.fv(K, data.frame(tcom=Ktrans), "bold(C)~hat(%s)[trans](r)",
                 nzpaste("model compensator of",
                         algo,
                         "translation-corrected %s"),
                 "trans")
  }
  if(opt$ripl) {
    # Ripley isotropic correction
    edgewt <- edge.Ripley(UI, matrix(dIJ, ncol=1))
    wh   <- whist(dIJ[okI],     breaks$val, (edgewt * wcIJ)[okI])
    whDD <- whist(dIJ[DD & okI], breaks$val, edgewt[DD & okI])    
#    Kiso <- cumsum(whDD)/(lambda2.used * area.used)
    Kiso <- cumsum(whDD)/(lambda2 * area.used)
    Kiso[r >= rmax] <- NA
    K <- bind.fv(K, data.frame(iso=Kiso), "hat(%s)[iso](r)",
                 nzpaste(algo,
                         "isotropic-corrected nonparametric estimate of %s"),
                 "iso")
#    lambda2U <- (npts.used + 1) * npts.used/(area.used^2)
    lambda2U <- (npts + 1) * npts/(areaW^2)    
    Kiso <- cumsum(wh)/(lambda2U * area.used)
    Kiso[r >= rmax] <- NA
    K <- bind.fv(K, data.frame(icom=Kiso), "bold(C)~hat(%s)[iso](r)",
                 nzpaste("model compensator of",
                         algo, "isotropic-corrected %s"),
                 "iso")
    #
    if(compute.var) {
      savedotnames <- fvnames(K, ".")
      # compute contribution to compensator from each quadrature point
      dOK <- dIJ[okI]
      eOK <- edgewt[okI]
      iOK <- I[okI]
      denom <- lambda2U * area.used
      variso <- varsumiso <- 0 * Kiso
      for(i in sort(unique(iOK))) {
        relevant <- (iOK == i)
        tincrem <- whist(dOK[relevant], breaks$val, eOK[relevant])
        localterm <- cumsum(tincrem)/denom
        variso <- variso + wc[i] * localterm^2
        if(Z[i])
          varsumiso <- varsumiso + localterm^2
      }
      sdiso <- sqrt(variso)
      K <- bind.fv(K, data.frame(ivar=variso,
                                 isd =sdiso,
                                 ihi = 2*sdiso,
                                 ilo = -2*sdiso,
                                 ivarsum=varsumiso),
                   c("bold(C)^2~hat(%s)[iso](r)",
                     "sqrt(bold(C)^2~hat(%s)[iso](r))",
                     "bold(R)~hat(%s)[hi](r)",
                     "bold(R)~hat(%s)[lo](r)",
                     "hat(C)^2~hat(%s)[iso](r)"),
                   c("Poincare variance of isotropic-corrected %s",
                     "sqrt(Poincare variance)  of isotropic-corrected %s",
                     "upper critical band for isotropic-corrected %s",
                     "lower critical band for isotropic-corrected %s",
                     "data estimate of Poincare variance of %s"),
                   "iso")
      # fvnames(K, ".") <- c(savedotnames, "isd")
      fvnames(K, ".") <- savedotnames
    }
  }

  # default is to display all corrections
  formula(K) <- . ~ r
  unitname(K) <- unitname(X)
  # secret tag used by 'Kres'
  attr(K, "maker") <- "Kcom"
  return(K)
}

# `reweighted' translation edge correction
edge.Trans.modif <- function(X, Y=X, WX=X$window, WY=Y$window,
                             exact=FALSE, paired=FALSE,
                             trim=spatstat.options("maxedgewt")) {

  # computes edge correction factor
  #  f = area(WY)/area(intersect.owin(WY, shift(WX, X[i] - Y[j])))
  
  X <- as.ppp(X, WX)

  W <- X$window
  x <- X$x
  y <- X$y

  Y <- as.ppp(Y, WY)
  xx <- Y$x
  yy <- Y$y

  nX <- npoints(X)
  nY <- npoints(Y)
  if(paired && (nX != nY))
    stop("X and Y should have equal length when paired=TRUE")
  
  # For irregular polygons, exact evaluation is very slow;
  # so use pixel approximation, unless exact=TRUE
  if(!exact) {
    if(WX$type == "polygonal")
      WX <- as.mask(WX)
    if(WY$type == "polygonal")
      WY <- as.mask(WX)
  }

  typeX <- WX$type
  typeY <- WY$type

  if(typeX == "rectangle" && typeY == "rectangle") {
    # Fast code for this case
    if(!paired) {
      DX <- abs(outer(x,xx,"-"))
      DY <- abs(outer(y,yy,"-"))
    } else {
      DX <- abs(xx - x)
      DY <- abs(yy - y)
    }
    A <- WX$xrange
    B <- WX$yrange
    a <- WY$xrange 
    b <- WY$yrange
    # compute width and height of intersection
    wide  <- pmin.int(a[2], A[2]+DX) - pmax(a[1], A[1]+DX)
    high  <- pmin.int(b[2], B[2]+DY) - pmax(b[1], B[1]+DY)
    # edge correction weight
    weight <- diff(a) * diff(b) / (wide * high)
    if(!paired)
      weight <- matrix(weight, nrow=X$n, ncol=Y$n)
  } else if(typeX %in% c("rectangle", "polygonal")
            && typeY %in% c("rectangle", "polygonal")) {
    # This code is SLOW
    WX <- as.polygonal(WX)
    WY <- as.polygonal(WY)
    a <- area(W)
    if(!paired) {
      weight <- matrix(, nrow=nX, ncol=nY)
      if(nX > 0 && nY > 0) {
        for(i in seq_len(nX)) {
          X.i <- c(x[i], y[i])
          for(j in seq_len(nY)) {
            shiftvector <- X.i - c(xx[j],yy[j])
            WXshift <- shift(WX, shiftvector)
            b <- overlap.owin(WY, WXshift)
            weight[i,j] <- a/b
          }
        }
      }
    } else {
      nX <- npoints(X)
      weight <- numeric(nX)
      if(nX > 0) {
        for(i in seq_len(nX)) {
          shiftvector <- c(x[i],y[i]) - c(xx[i],yy[i])
          WXshift <- shift(WX, shiftvector)
          b <- overlap.owin(WY, WXshift)
          weight[i] <- a/b
        }
      }
    }
  } else {
    WX <- as.mask(WX)
    WY <- as.mask(WY)
    # make difference vectors
    if(!paired) {
      DX <- outer(x,xx,"-")
      DY <- outer(y,yy,"-")
    } else {
      DX <- x - xx
      DY <- y - yy
    }
    # compute set cross-covariance
    g <- setcov(WY,WX)
    # evaluate set cross-covariance at these vectors
    gvalues <- lookup.im(g, as.vector(DX), as.vector(DY),
                         naok=TRUE, strict=FALSE)
    weight <- area(WY)/gvalues
  }

  # clip high values
  if(length(weight) > 0)
    weight <- pmin.int(weight, trim)
  if(!paired) 
    weight <- matrix(weight, nrow=X$n, ncol=Y$n)
  return(weight)
}

Kcom
})
