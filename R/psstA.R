#
#	psstA.R
#
#	Pseudoscore residual for unnormalised F (area-interaction)
#
#	$Revision: 1.7 $	$Date: 2014/11/11 02:31:44 $
#
################################################################################
#

psstA <- function(object, r=NULL, breaks=NULL, ...,
                  model=NULL,
                  trend=~1, interaction=Poisson(),
                  rbord=reach(interaction), ppmcorrection="border",
                  correction="all",
                  truecoef=NULL, hi.res=NULL,
                  nr=spatstat.options("psstA.nr"),
                  ngrid=spatstat.options("psstA.ngrid")) {
  if(is.ppm(object))
    fit <- object
  else if(is.ppp(object) || is.quad(object)) {
    # convert to quadscheme
    if(is.ppp(object))
      object <- quadscheme(object, ...)
    # fit model
    if(!is.null(model))
      fit <- update(model, Q=object, forcefit=TRUE)
    else if(ppmcorrection == "border")
      fit <- ppm(object,
                 trend=trend, interaction=interaction,
                 rbord=rbord, forcefit=TRUE)
    else
      fit <- ppm(object,
                 trend=trend, interaction=interaction,
                 correction=ppmcorrection, forcefit=TRUE)
  } else 
    stop("object should be a fitted point process model or a point pattern")

  rfixed <- !is.null(r) || !is.null(breaks)
  
  # Extract data and quadrature points
  Q <- quad.ppm(fit, drop=FALSE)
  X <- data.ppm(fit)
  U <- union.quad(Q)
  Z <- is.data(Q) # indicator data/dummy
#  E <- equalsfun.quad(Q)
#  WQ <- w.quad(Q)  # quadrature weights

  # integrals will be restricted to quadrature points
  # that were actually used in the fit
#  USED <- getglmsubset(fit)
  if(fit$correction == "border") {
    rbord <- fit$rbord
    b <- bdist.points(U)
    USED <- (b > rbord)
    bX <- bdist.points(X)
    USEDX <- (bX > rbord)
  } else {
    USED <- rep.int(TRUE, U$n)
    USEDX <- rep.int(TRUE, X$n)
  }
  
  # basic statistics
  Win <- Window(X)
  npts <- npoints(X)
  areaW <- area(Win)
  lambda <- npts/areaW

  #  determine breakpoints for r values
  rmaxdefault <- rmax.rule("F", Win, lambda)
  if(rfixed) 
    breaks <- handle.r.b.args(r, breaks, Win, rmaxdefault=rmaxdefault)
  else {
    # create fairly coarse 'r' values
    r <- seq(0, rmaxdefault, length=nr)
    breaks <- breakpts.from.r(r)
  }
  rvals <- breaks$r
  rmax  <- breaks$max
  
  # residuals
  res <- residuals(fit, type="raw", drop=FALSE,
                    new.coef=truecoef, quad=hi.res)
  # 
  rescts <- with(res, "continuous")
  # absolute weight for continuous integrals
  wc   <- -rescts

  # initialise fv object
  df <- data.frame(r=rvals, theo=0)
  desc <- c("distance argument r", "value 0 corresponding to perfect fit")
  ans <- fv(df, "r", substitute(bold(R)~Delta~V[A](r), NULL),
            "theo", . ~ r,
            alim=c(0, rmax), c("r","%s[theo](r)"), desc,
            fname="bold(R)~Delta~V[A]")

  #
  # for efficiency, compute the largest value of distance transform
  Dmax <- 0
  for(i in 1:npts) {
    Di <- distmap(X[-i])
    Dimax <- summary(Di)$max
    Dmax <- max(Dmax, Dimax)
  }
  Rmax <- min(max(rvals), Dmax * 1.1)
  nontrivial <- (rvals <= Rmax)
  trivialzeroes <- numeric(sum(!nontrivial))
  
  # pseudosum
  Ax <- areaLoss.grid(X, rvals[nontrivial], subset=USEDX, ngrid=ngrid)
  C1 <- apply(Ax, 2, sum)
  C1 <- c(C1, trivialzeroes)
  # pseudocompensator
  OK <- USED & !Z
  Au <- areaGain.grid(U[OK], X, rvals[nontrivial], W=Win, ngrid=ngrid)
  lamu <- matrix(wc[OK], nrow=nrow(Au), ncol=ncol(Au))
  C2 <- apply(lamu * Au, 2, sum)
  C2 <- c(C2, trivialzeroes)
  # pseudoscore residual
  Ctot <- C1 - C2
  # tack on
  ans <- bind.fv(ans,
                 data.frame(dat=C1,
                            com=C2,
                            res=Ctot),
                 c("Sigma~Delta~V[A](r)", "bold(C)~Delta~V[A](r)", "%s(r)"),
                 c("data pseudosum (contribution to %s)",
                   "model pseudocompensator (contribution to %s)",
                   "pseudoscore residual %s"),
               "res")
  #
  # pseudovariance
  #        (skipped if called by envelope() etc)
  #
  if(correction == "all") {
    lamX <- matrix(wc[USED & Z], nrow=nrow(Ax), ncol=ncol(Ax))
    Var <- apply(lamu * Au^2, 2, sum) + apply(lamX * Ax^2, 2, sum)
    Var <- c(Var, trivialzeroes)
    # two-sigma limits
    TwoSig <- 2 * sqrt(Var)
    # tack on
    ans <- bind.fv(ans,
                   data.frame(var=Var,
                              up=TwoSig,
                              lo=-TwoSig),
                 c("bold(C)^2~Delta~V[A](r)",
                   "%s[up](r)", "%s[lo](r)"),
                 c("pseudovariance of %s",
                   "upper 2sigma critical limit for %s",
                   "lower 2sigma critical limit for %s"),
               "res")
    fvnames(ans, ".") <- c("res", "up", "lo", "theo")
  }
  unitname(ans) <- unitname(fit)
  # 
  return(ans)
}
