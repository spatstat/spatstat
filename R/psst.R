#
#	psst.R
#
#	Computes the GNZ contrast of delta-f for any function f
#
#	$Revision: 1.5 $	$Date: 2013/04/25 06:37:43 $
#
################################################################################
#

psst <- function(object, fun, r=NULL, breaks=NULL, ...,
                 trend=~1, interaction=Poisson(),
                 rbord=reach(interaction),
                 truecoef=NULL, hi.res=NULL,
                 funargs=list(correction="best"),
                 verbose=TRUE) {
  if(inherits(object, "ppm")) 
    fit <- object
  else if(inherits(object, "ppp"))
    fit <- ppm(quadscheme(object, ...),
               trend=trend, interaction=interaction, rbord=rbord)
  else if(inherits(object, "quad")) 
    fit <- ppm(object, 
               trend=trend, interaction=interaction, rbord=rbord)
  else 
    stop("object should be a fitted point process model or a point pattern")

  rfixed <- !is.null(r) || !is.null(breaks)
  
  # Extract data and quadrature points
  Q <- quad.ppm(fit, drop=FALSE)
  X <- data.ppm(fit)
  U <- union.quad(Q)
  Z <- is.data(Q) # indicator data/dummy
  E <- equalsfun.quad(Q)
  WQ <- w.quad(Q)  # quadrature weights

  # integrals will be restricted to quadrature points
  # that were actually used in the fit
#  USED <- getglmsubset(fit)
  if(fit$correction == "border") {
    rbord <- fit$rbord
    b <- bdist.points(U)
    USED <- (b > rbord)
  } else USED <- rep.int(TRUE, U$n)
  
  # basic statistics
  Win <- X$window
  npoints <- X$n
  area <- area.owin(Win)
  lambda <- npoints/area

  # adjustments to account for restricted domain of pseudolikelihood
  if(any(!USED)) {
    XUSED <- USED[Z]
    npoints.used <- sum(Z & USED)
    area.used <- sum(WQ[USED])
    lambda.used <- npoints.used/area.used
  } else {
    XUSED <- rep.int(TRUE, npoints)
    npoints.used <- npoints
    area.used <- area
    lambda.used <- lambda
  }
  
  #  determine breakpoints for r values
  rmaxdefault <- rmax.rule("G", Win, lambda)
  breaks <- handle.r.b.args(r, breaks, Win, rmaxdefault=rmaxdefault)
  rvals <- breaks$r
  rmax  <- breaks$max
  
  # residuals
  resid <- residuals(fit, type="raw",drop=FALSE,
                    coefs=truecoef, quad=hi.res)
  rescts <- with(resid, "continuous")
  # absolute weight for continuous integrals
  wc   <- -rescts

  # initialise fv object
  df <- data.frame(r=rvals, theo=0)
  desc <- c("distance argument r", "value 0 corresponding to perfect fit")
  ans <- fv(df, "r", substitute(bold(R)~Delta~S(r), NULL),
            "theo", . ~ r,
            alim=c(0, rmax), c("r","%s[theo](r)"), desc,
            fname="bold(R)~Delta~S")

  # evaluate fun(X) for data
  fX <- do.call(fun, append(list(X, r=rvals), funargs))
  fXunits <- unitname(fX)
  # Extract 'best' estimate only
  fX <- with(fX, .y)
  zero <- numeric(length(fX))
  # sum over all quadrature points
  iused <- seq(U$n)[USED]
  nused <- length(iused)
  if(verbose) cat(paste("\nProcessing", nused, "quadrature points..."))
  # running sums & integrals
  sumX <- zero
  integ <- integ2 <- zero
  # template for X \cup {u}
  uX <- superimpose(U[1], X, W=Win, check=FALSE)
  Ux <- U$x
  Uy <- U$y
  # 
  for(j in seq(nused)) {
    i <- iused[j]
    wi <- wc[i]
    if(Z[i]) {
      # data point
      fXi <- do.call(fun, append(list(X[-i], r=rvals), funargs))
      fXi <- with(fXi, .y)
      deltaf <- fX - fXi
      sumX <- sumX + deltaf
    } else {
      # dummy point
      uX$x[1] <- Ux[i]
      uX$y[1] <- Uy[i]
      fuX <- do.call(fun, append(list(uX, r=rvals), funargs))
      fuX <- with(fuX, .y)
      deltaf <- fuX - fX
    }
    integ <- integ + wi * deltaf
    integ2 <- integ2 + wi * deltaf^2
    # 
    if(j %% 500 == 0) {
      cat("[garbage ")
      gc()
      cat("collected]")
    }
    if(verbose) progressreport(j, nused)
  }

  sdv <- sqrt(integ2)
  res <- sumX - integ
  ans <- bind.fv(ans,
                 data.frame(dat=sumX,
                            com=integ,
                            var=integ2,
                            sd=sdv,
                            hi=2*sdv,
                            lo=-2*sdv,
                            res=res,
                            stdres=res/sdv),
                 c("Sigma~Delta~S(r)",
                   "bold(C)~Delta~S(r)",
                   "bold(C)^2~Delta~S(r)",
                   "sqrt(bold(C)^2~Delta~S(r))",
                   "%s[hi](r)",
                   "%s[lo](r)",
                   "bold(R)~Delta~S(r)",
                   "bold(T)~Delta~S(r)"),
               c("data pseudosum (contribution to %s)",
                 "model compensator (contribution to %s)",
                 "pseudovariance of %s",
                 "sqrt(pseudovariance) of %s",
                 "upper 2 sigma critical band for %s",
                 "lower 2 sigma critical band for %s",
                 "pseudoresidual function %s",
                 "standardised pseudoresidual function %s"),
               "res")

  fvnames(ans,".") <- c("res", "hi", "lo", "theo")
  unitname(ans) <- fXunits
  # 
  return(ans)
}

npfun <- function(X, ..., r) {
  npts <- npoints(X)
  # initialise fv object
  df <- data.frame(r=r, theo=0, npoint=npts)
  desc <- c("distance argument r",
            "value 0",
            "value equal to number of points")
  ans <- fv(df, "r", substitute(npoints(r), NULL),
            "npoint", . ~ r,
            alim=c(0, max(r)), c("r","%s[theo](r)", "%s[obs](r)"),
            desc, fname="npoints")
  unitname(ans) <- unitname(X)
  return(ans)
}

nndcumfun <- function(X, ..., r) {
  nn <- nndist(X)
  bk <- breakpts.from.r(r)
#  nn <- nn[nn <= bdist.points(X)]
  h <- whist(nn, bk$val)
  # initialise fv object
  df <- data.frame(r=r, theo=0, obs=h)
  desc <- c("distance argument r",
            "value 0",
            "observed count")
  ans <- fv(df, "r", substitute(nndcount(r), NULL),
            "obs", . ~ r,
            alim=c(0, max(r)), c("r","%s[theo](r)", "%s[obs](r)"),
            desc, fname="nndcount")
  unitname(ans) <- unitname(X)
  return(ans)
}

  
