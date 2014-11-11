#
#	psstG.R
#
#	Pseudoscore residual for unnormalised G (saturation process)
#
#	$Revision: 1.7 $	$Date: 2014/11/11 02:33:16 $
#
################################################################################
#

psstG <- function(object, r=NULL, breaks=NULL, ...,
                  model=NULL, 
                  trend=~1, interaction=Poisson(),
                  rbord=reach(interaction),
                  truecoef=NULL, hi.res=NULL) {
  if(inherits(object, "ppm")) 
    fit <- object
  else if(inherits(object, "ppp") || inherits(object, "quad")) {
    # convert to quadscheme
    if(inherits(object, "ppp"))
      object <- quadscheme(object, ...)
    # fit model
    if(!is.null(model))
      fit <- update(model, Q=object, forcefit=TRUE)
    else 
      fit <- ppm(object,
                 trend=trend, interaction=interaction,
                 rbord=rbord, forcefit=TRUE)
  } else 
    stop("object should be a fitted point process model or a point pattern")

#  rfixed <- !is.null(r) || !is.null(breaks)
  
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
  Win <- Window(X)
  npts <- npoints(X)
  areaW <- area(Win)
  lambda <- npts/areaW

  # adjustments to account for restricted domain of pseudolikelihood
#  if(any(!USED)) {
#    npts.used <- sum(Z & USED)
#    area.used <- sum(WQ[USED])
#    lambda.used <- npts.used/area.used
#  } else {
#    npts.used <- npts
#    area.used <- areaW
#    lambda.used <- lambda
#  }
  
  #  determine breakpoints for r values
  rmaxdefault <- rmax.rule("G", Win, lambda)
  breaks <- handle.r.b.args(r, breaks, Win, rmaxdefault=rmaxdefault)
  rvals <- breaks$r
  rmax  <- breaks$max
  
  # residuals
  res <- residuals(fit, type="raw",drop=FALSE,
                    new.coef=truecoef, quad=hi.res)
  resval <- with(res, "increment")
  rescts <- with(res, "continuous")
  # absolute weight for continuous integrals
  wc   <- -rescts

  # initialise fv object
  df <- data.frame(r=rvals, theo=0)
  desc <- c("distance argument r", "value 0 corresponding to perfect fit")
  ans <- fv(df, "r", substitute(bold(R)~Delta~V[S](r), NULL),
            "theo", . ~ r,
            alim=c(0, rmax), c("r","%s[theo](r)"), desc,
            fname="bold(R)~Delta~V[S]")

  # First phase: .................................................
  # nearest neighbours (quadrature point to data point)
  nn <- nncross(U, X, seq(U$n), seq(X$n)) # excludes identical pairs
  dIJ <- nn$dist
  I <- seq(U$n)
  J <- nn$which
  DD <- (I <= X$n)  # TRUE for data points
  wcIJ <- wc
  okI <- USED[I]

  # histogram of nndist for data points only (without edge correction)
  Bsum <- cumsum(whist(dIJ[DD & okI], breaks$val))
  # weighted histogram of nncross (without edge correction)
  Bint <- cumsum(whist(dIJ[okI], breaks$val, wcIJ[okI]))
  # residual
  Bres <- Bsum - Bint
  # tack on 
  
  ans <- bind.fv(ans,
                 data.frame(dat1=Bsum,
                            com1=Bint,
                            res1=Bres),
                 c("%s[dat1](r)",
                   "%s[com1](r)",
                   "%s[res1](r)"),
                 c("phase 1 pseudosum (contribution to %s)",
                   "phase 1 pseudocompensator (contribution to %s)",
                   "phase 1 pseudoresidual (contribution to %s)"))
  
  # Second phase: ................................................
  # close pairs (quadrature point to data point)
  close <- crosspairs(U, X, rmax)
  dIJ <- close$d
  I   <- close$i
  J   <- close$j
#  UI <- U[I]
#  XJ <- X[J]
  EIJ <- E(I, J) # TRUE if points are identical, U[I[k]] == X[J[k]] 
  ZI <- Z[I]     # TRUE if U[I[k]] is a data point
  DD <- ZI & !EIJ  # TRUE for pairs of distinct data points only
#  nDD <- sum(DD)
  okI <- USED[I]
  
  # residual weights
#  wIJ <- ifelseXY(EIJ, rescts[I], resval[I])
  # absolute weight for continuous integrals
  wc   <- -rescts
  wcIJ <- -rescts[I]
  
  # nearest and second-nearest neighbour distances in X
  nn1 <- nndist(X)
  nn2 <- nndist(X, k=2)
  nn1J <- nn1[J]
  nn2J <- nn2[J]
  
  # weird use of the reduced sample estimator
  # data sum:
  RSX <- Kount(dIJ[DD & okI], nn2J[DD & okI], nn2J[ZI & okI], breaks)
  Csum <- RSX$numerator
  # integral:
  if(spatstat.options("psstG.remove.zeroes"))
    okE <- okI & !EIJ
  else
    okE <- okI
  RSD <- Kwtsum(dIJ[okE], nn1J[okE], wcIJ[okE],
                  nn1, rep.int(1, length(nn1)), breaks)
  Cint <- RSD$numerator
  #
  Cres <- Bres + Csum - Cint
  # tack on 
  ans <- bind.fv(ans,
                 data.frame(dat2=Csum,
                            com2=Cint,
                            res2=Cres,
                            dat=Bsum+Csum,
                            com=Bint+Cint,
                            res=Bres+Cres),
                 c("%s[dat2](r)",
                   "%s[com2](r)",
                   "%s[res2](r)",
                   "Sigma~Delta~V[S](r)",
                   "bold(C)~Delta~V[S](r)",
                   "bold(R)~Delta~V[S](r)"),
                 c("phase 2 pseudosum (contribution to %s)",
                   "phase 2 pseudocompensator (contribution to %s)",
                   "phase 2 pseudoresidual (contribution to %s)",
                   "pseudosum (contribution to %s)",
                   "pseudocompensator (contribution to %s)",
                   "pseudoresidual function %s"),
                 "res")
  # restrict choice of curves in default plot
  fvnames(ans, ".") <- c("dat", "com", "res", "theo")
  # 
  return(ans)
}
