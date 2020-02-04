#'
#'    tests/deltasuffstat.R
#' 
#'    Explicit tests of 'deltasuffstat'
#' 
#' $Revision: 1.2 $ $Date: 2018/03/28 09:04:06 $

require(spatstat)
local({
  
  disagree <- function(x, y, tol=1e-7) { !is.null(x) && !is.null(y) && max(abs(x-y)) > tol }

  flydelta <- function(model, modelname="") {
    ## Check execution of different algorithms for 'deltasuffstat' 
    dSS <- deltasuffstat(model, sparseOK=TRUE)
    dBS <- deltasuffstat(model, sparseOK=TRUE,  use.special=FALSE, force=TRUE)
    dBF <- deltasuffstat(model, sparseOK=FALSE, use.special=FALSE, force=TRUE)
    ## Compare results
    if(disagree(dBS, dSS))
      stop(paste(modelname, "model: Brute force algorithm disagrees with special algorithm"))
    if(disagree(dBF, dBS))
      stop(paste(modelname, "model: Sparse and full versions of brute force algorithm disagree"))
    return(invisible(NULL))
  }

  modelS <- ppm(cells ~ x, Strauss(0.13), nd=10)
  flydelta(modelS, "Strauss")

  antsub <- ants[c(FALSE,TRUE,FALSE)]
  rmat <- matrix(c(130, 90, 90, 60), 2, 2)
  
  modelM <- ppm(antsub ~ 1, MultiStrauss(rmat), nd=16)
  flydelta(modelM, "MultiStrauss")
                
  modelA <- ppm(antsub ~ 1, HierStrauss(rmat, archy=c(2,1)), nd=16)
  flydelta(modelA, "HierStrauss")
})
#'
#'  tests/density.R
#'
#'  Test behaviour of density() methods,
#'                    relrisk(), Smooth()
#'                    and inhomogeneous summary functions
#'                    and idw, adaptive.density, intensity
#'
#'  $Revision: 1.52 $  $Date: 2020/02/04 01:54:56 $
#'

require(spatstat)

local({

  # test all cases of density.ppp and densityfun.ppp
  
  tryit <- function(..., do.fun=TRUE) {
    Z <- density(cells, ..., at="pixels")
    Z <- density(cells, ..., at="points")
    if(do.fun) {
      f <- densityfun(cells, ...)
      U <- f(0.1, 0.3)
    }
    return(invisible(NULL))
  }

  tryit(0.05)
  tryit(0.05, diggle=TRUE)
  tryit(0.05, se=TRUE)
  tryit(0.05, weights=expression(x))

  tryit(0.07, kernel="epa")
  tryit(0.07, kernel="quartic")
  tryit(0.07, kernel="disc")

  tryit(0.07, kernel="epa", weights=expression(x))

  tryit(sigma=Inf)
  tryit(sigma=Inf, weights=expression(x))
  
  V <- diag(c(0.05^2, 0.07^2))
  tryit(varcov=V)
  tryit(varcov=V, diggle=TRUE)
  tryit(varcov=V, weights=expression(x))
  tryit(varcov=V, weights=expression(x), diggle=TRUE)
  Z <- distmap(runifpoint(5, Window(cells)))
  tryit(0.05, weights=Z)
  tryit(0.05, weights=Z, diggle=TRUE)

  trymost <- function(...) tryit(..., do.fun=FALSE) 
  wdf <- data.frame(a=1:42,b=42:1)
  trymost(0.05, weights=wdf)
  trymost(0.05, weights=wdf, diggle=TRUE)
  trymost(sigma=Inf, weights=wdf)
  trymost(varcov=V, weights=wdf)
  trymost(varcov=V, weights=expression(cbind(x,y)))
  
  ## run C algorithm 'denspt'
  opa <- spatstat.options(densityC=TRUE, densityTransform=FALSE)
  tryit(varcov=V)
  tryit(varcov=V, weights=expression(x))
  trymost(varcov=V, weights=wdf)
  spatstat.options(opa)
  
  crossit <- function(..., sigma=NULL) {
    U <- runifpoint(20, Window(cells))
    a <- densitycrossEngine(cells, U, ..., sigma=sigma)
    a <- densitycrossEngine(cells, U, ..., sigma=sigma, diggle=TRUE)
    invisible(NULL)
  }
  crossit(varcov=V, weights=cells$x)
  crossit(varcov=V, weights=wdf)
  crossit(sigma=0.1, weights=wdf)
  crossit(sigma=0.1, kernel="epa", weights=wdf)

  crossit(sigma=Inf)
  
  # apply different discretisation rules
  Z <- density(cells, 0.05, fractional=TRUE)
  Z <- density(cells, 0.05, preserve=TRUE)
  Z <- density(cells, 0.05, fractional=TRUE, preserve=TRUE)
        
  ## compare results with different algorithms
  crosscheque <- function(expr) {
    e <- as.expression(substitute(expr))
    ename <- sQuote(deparse(substitute(expr)))
    ## interpreted R
    opa <- spatstat.options(densityC=FALSE, densityTransform=FALSE)
    val.interpreted <- eval(e)
    ## established C algorithm 'denspt'
    spatstat.options(densityC=TRUE, densityTransform=FALSE)
    val.C <- eval(e)
    ## new C algorithm 'Gdenspt' using transformed coordinates
    spatstat.options(densityC=TRUE, densityTransform=TRUE)
    val.Transform <- eval(e)
    spatstat.options(opa)
    if(max(abs(val.interpreted - val.C)) > 0.001)
      stop(paste("Numerical discrepancy between R and C algorithms in",
                 ename))
    if(max(abs(val.C - val.Transform)) > 0.001)
      stop(paste("Numerical discrepancy between C algorithms",
                 "using transformed and untransformed coordinates in",
                 ename))
    invisible(NULL)
  }

  ## execute & compare results of density(at="points") with different algorithms
  wdfr <- cbind(1:npoints(redwood), 2)
  crosscheque(density(redwood, at="points", sigma=0.13, edge=FALSE))
  crosscheque(density(redwood, at="points", sigma=0.13, edge=FALSE,
                      weights=wdfr[,1]))
  crosscheque(density(redwood, at="points", sigma=0.13, edge=FALSE,
                      weights=wdfr))

  ## correctness of non-Gaussian kernel calculation
  leavein <- function(ker, maxd=0.025) {
    ZI <- density(redwood, 0.12, kernel=ker, edge=FALSE,
                  dimyx=256)[redwood]
    ZP <- density(redwood, 0.12, kernel=ker, edge=FALSE,
                  at="points", leaveoneout=FALSE)
    discrep <- max(abs(ZP - ZI))/npoints(redwood)
    if(discrep > maxd) 
      stop(paste("Discrepancy",
                 signif(discrep, 3),
                 "in calculation for", ker, "kernel"))
    return(invisible(NULL))
  }
  leavein("epanechnikov", 0.015)
  leavein("quartic",      0.010)
  leavein("disc",         0.100) 

  ## bandwidth selection code blocks
  sigvec <- 0.01 * 2:15
  sigran <- range(sigvec)
  bw.ppl(redwood, sigma=sigvec)
  bw.ppl(redwood, srange=sigran, ns=5)
  bw.CvL(redwood, sigma=sigvec)
  bw.CvL(redwood, srange=sigran, ns=5)

  ## adaptive bandwidth
  a <- bw.abram(redwood)
  a <- bw.abram(redwood, pilot=density(redwood, 0.2))
  a <- bw.abram(redwood, smoother="densityVoronoi", at="pixels")
  
  ## Kinhom 
  lam <- density(redwood)
  K <- Kinhom(redwood, lam)
  
  lamX <- density(redwood, at="points")
  KX <- Kinhom(redwood, lamX)

  ## test all code cases of new 'relrisk.ppp' algorithm
  pants <- function(..., X=ants, sigma=100, se=TRUE) {
    a <- relrisk(X, sigma=sigma, se=se, ...)
    return(TRUE)
  }
  pants()
  pants(diggle=TRUE)
  pants(edge=FALSE)
  pants(diggle=TRUE, at="points")
  pants(edge=FALSE, at="points")
  pants(casecontrol=FALSE)
  pants(relative=TRUE)
  pants(casecontrol=FALSE, relative=TRUE)
  pants(at="points")
  pants(casecontrol=FALSE,at="points")
  pants(relative=TRUE,at="points")
  pants(casecontrol=FALSE, relative=TRUE,at="points")
  pants(relative=TRUE, control="Cataglyphis", case="Messor")
  pants(relative=TRUE, control="Cataglyphis", case="Messor", at="points")
  pants(casecontrol=FALSE, case="Messor", se=FALSE)
  pants(case=2, at="pixels", relative=TRUE)
  pants(case=2, at="points", relative=TRUE)
  pants(case=2, at="pixels", relative=FALSE)
  pants(case=2, at="points", relative=FALSE)
  pants(sigma=Inf)
  pants(sigma=NULL, varcov=diag(c(100,100)^2))
  
  ## more than 2 types
  pants(X=sporophores)
  pants(X=sporophores, sigma=20, at="points")
  pants(X=sporophores, sigma=20, relative=TRUE, at="points")
  pants(X=sporophores, sigma=20, at="pixels", se=FALSE)
  pants(X=sporophores, sigma=20, relative=TRUE, at="pixels", se=FALSE)
  bw.relrisk(sporophores, method="leastsquares")
  bw.relrisk(sporophores, method="weightedleastsquares")
  
  ## likewise 'relrisk.ppm'
  fit <- ppm(ants ~ x)
  rants <- function(..., model=fit) {
    a <- relrisk(model, sigma=100, se=TRUE, ...)
    return(TRUE)
  }
  rants()
  rants(diggle=TRUE)
  rants(edge=FALSE)
  rants(diggle=TRUE, at="points")
  rants(edge=FALSE, at="points")
  rants(casecontrol=FALSE)
  rants(relative=TRUE)
  rants(casecontrol=FALSE, relative=TRUE)
  rants(at="points")
  rants(casecontrol=FALSE,at="points")
  rants(relative=TRUE,at="points")
  rants(casecontrol=FALSE, relative=TRUE,at="points")
  rants(relative=TRUE, control="Cataglyphis", case="Messor")
  rants(relative=TRUE, control="Cataglyphis", case="Messor", at="points")

  ## more than 2 types
  fut <- ppm(sporophores ~ x)
  rants(model=fut)
  rants(model=fut, at="points")
  rants(model=fut, relative=TRUE, at="points")
  
  ## execute Smooth.ppp and Smoothfun.ppp in all cases
  stroke <- function(..., Y = longleaf) {
    Z <- Smooth(Y, ..., at="pixels")
    Z <- Smooth(Y, ..., at="points", leaveoneout=TRUE)
    Z <- Smooth(Y, ..., at="points", leaveoneout=FALSE)
    f <- Smoothfun(Y, ...)
    f(120, 80)
    f(Y[1:2])
    f(Y[FALSE])
    U <- as.im(f)
    return(invisible(NULL))
  }
  stroke()
  stroke(5, diggle=TRUE)
  stroke(5, geometric=TRUE)
  stroke(1e-6) # generates warning about small bandwidth
  stroke(5, weights=runif(npoints(longleaf)))
  stroke(5, weights=expression(x))
  stroke(5, kernel="epa")
  stroke(varcov=diag(c(25, 36)))
  stroke(varcov=diag(c(25, 36)), weights=runif(npoints(longleaf)))
  stroke(5, Y=longleaf %mark% 1)
  stroke(5, Y=cut(longleaf,breaks=3))
  Z <- as.im(function(x,y){abs(x)+1}, Window(longleaf))
  stroke(5, weights=Z)
  stroke(5, weights=Z, geometric=TRUE)

  stroke(sigma=Inf)
  
  markmean(longleaf, 9)
  
  strike <- function(..., Y=finpines) {
    Z <- Smooth(Y, ..., at="pixels")
    Z <- Smooth(Y, ..., at="points", leaveoneout=TRUE)
    Z <- Smooth(Y, ..., at="points", leaveoneout=FALSE)
    f <- Smoothfun(Y, ...)
    f(4, 1)
    f(Y[1:2])
    f(Y[FALSE])
    U <- as.im(f)
    return(invisible(NULL))
  }
  strike()
  strike(sigma=1.5, kernel="epa")
  strike(varcov=diag(c(1.2, 2.1)))
  strike(sigma=1e-6)
  strike(sigma=1e-6, kernel="epa")
  strike(sigma=Inf)
  strike(1.5, weights=runif(npoints(finpines)))
  strike(1.5, weights=expression(y))
  strike(1.5, geometric=TRUE)
  strike(1.5, Y=finpines[FALSE])
  flatfin <- finpines %mark% data.frame(a=rep(1, npoints(finpines)), b=2)
  strike(1.5, Y=flatfin)
  strike(1.5, Y=flatfin, geometric=TRUE)

  opx <- spatstat.options(densityTransform=FALSE)
  stroke(5, Y=longleaf[order(longleaf$x)], sorted=TRUE)
  strike(1.5, Y=finpines[order(finpines$x)], sorted=TRUE)
  spatstat.options(opx)

  ## detect special cases
  Smooth(longleaf[FALSE])
  Smooth(longleaf, minnndist(longleaf))
  Xconst <- cells %mark% 1
  Smooth(Xconst, 0.1)
  Smooth(Xconst, 0.1, at="points")
  Smooth(cells %mark% runif(42), sigma=Inf)
  Smooth(cells %mark% runif(42), sigma=Inf, at="points")
  Smooth(cells %mark% runif(42), sigma=Inf, at="points", leaveoneout=FALSE)
  Smooth(cut(longleaf, breaks=4))

  ## code not otherwise reached
  smoothpointsEngine(cells, values=rep(1, npoints(cells)), sigma=0.2)
  smoothpointsEngine(cells, values=runif(npoints(cells)), sigma=Inf)
  smoothpointsEngine(cells, values=runif(npoints(cells)), sigma=1e-16)
  
  ## validity of Smooth.ppp(at='points')
  Y <- longleaf %mark% runif(npoints(longleaf), min=41, max=43)
  Z <- Smooth(Y, 5, at="points", leaveoneout=TRUE)
  rZ <- range(Z)
  if(rZ[1] < 40 || rZ[2] > 44)
    stop("Implausible results from Smooth.ppp(at=points, leaveoneout=TRUE)")

  Z <- Smooth(Y, 5, at="points", leaveoneout=FALSE)
  rZ <- range(Z)
  if(rZ[1] < 40 || rZ[2] > 44)
    stop("Implausible results from Smooth.ppp(at=points, leaveoneout=FALSE)")

  ## compare Smooth.ppp results with different algorithms
  crosscheque(Smooth(longleaf, at="points", sigma=6))
  wt <- runif(npoints(longleaf))
  vc <- diag(c(25,36))
  crosscheque(Smooth(longleaf, at="points", sigma=6, weights=wt))
  crosscheque(Smooth(longleaf, at="points", varcov=vc))
  crosscheque(Smooth(longleaf, at="points", varcov=vc, weights=wt))

  ## drop-dimension coding errors
  X <- longleaf
  marks(X) <- cbind(marks(X), 1)
  Z <- Smooth(X, 5)

  ZZ <- bw.smoothppp(finpines, hmin=0.01, hmax=0.012, nh=2) # reshaping problem

  ## geometric-mean smoothing
  U <- Smooth(longleaf, 5, geometric=TRUE)
  UU <- Smooth(X, 5, geometric=TRUE)
  V <- Smooth(longleaf, 5, geometric=TRUE, at="points")
  VV <- Smooth(X, 5, geometric=TRUE, at="points")
})

reset.spatstat.options()

local({
  #' Kmeasure, second.moment.engine
  #' Expansion of window
  Zno  <- Kmeasure(redwood, sigma=0.2, expand=FALSE)
  Zyes <- Kmeasure(redwood, sigma=0.2, expand=TRUE)
  #' All code blocks
  sigmadouble <- rep(0.1, 2)
  diagmat <- diag(sigmadouble^2)
  generalmat <- matrix(c(1, 0.5, 0.5, 1)/100, 2, 2)
  Z <- Kmeasure(redwood, sigma=sigmadouble)
  Z <- Kmeasure(redwood, varcov=diagmat)
  Z <- Kmeasure(redwood, varcov=generalmat)
  A <- second.moment.calc(redwood, 0.1, what="all", debug=TRUE)
  B <- second.moment.calc(redwood, varcov=diagmat,    what="all")
  B <- second.moment.calc(redwood, varcov=diagmat,    what="all")
  D <- second.moment.calc(redwood, varcov=generalmat, what="all")
  PR <- pixellate(redwood)
  DRno  <- second.moment.calc(PR, 0.2, debug=TRUE, expand=FALSE,
                              npts=npoints(redwood), obswin=Window(redwood))
  DRyes <- second.moment.calc(PR, 0.2, debug=TRUE, expand=TRUE,
                              npts=npoints(redwood), obswin=Window(redwood))
  DR2 <- second.moment.calc(solist(PR, PR), 0.2, debug=TRUE, expand=TRUE,
                            npts=npoints(redwood), obswin=Window(redwood))
  Gmat <- generalmat * 100
  isoGauss <- function(x,y) {dnorm(x) * dnorm(y)}
  ee <- evaluate2Dkernel(isoGauss, runif(10), runif(10),
                         varcov=Gmat, scalekernel=TRUE)
  isoGaussIm <- as.im(isoGauss, square(c(-3,3)))
  gg <- evaluate2Dkernel(isoGaussIm, runif(10), runif(10),
                         varcov=Gmat, scalekernel=TRUE)
  ## experimental code
  op <- spatstat.options(developer=TRUE)
  DR <- density(redwood, 0.1)
  spatstat.options(op)
})

local({
  #' bandwidth selection
  op <- spatstat.options(n.bandwidth=8)
  bw.diggle(cells) 
  bw.diggle(cells, method="interpreted") # undocumented test
#  bw.relrisk(urkiola, hmax=20) is tested in man/bw.relrisk.Rd
  bw.relrisk(urkiola, hmax=20, method="leastsquares")
  bw.relrisk(urkiola, hmax=20, method="weightedleastsquares")
  ZX <- density(swedishpines, at="points")
  bw.pcf(swedishpines, lambda=ZX)
  bw.pcf(swedishpines, lambda=ZX,
         bias.correct=FALSE, simple=FALSE, cv.method="leastSQ")
  spatstat.options(op)
})

local({
  #' code in kernels.R
  kernames <- c("gaussian", "rectangular", "triangular",
                "epanechnikov", "biweight", "cosine", "optcosine")
  X <- rnorm(20)
  U <- runif(20)
  for(ker in kernames) {
    dX <- dkernel(X, ker)
    fX <- pkernel(X, ker)
    qU <- qkernel(U, ker)
    m0 <- kernel.moment(0, 0, ker)
    m1 <- kernel.moment(1, 0, ker)
    m2 <- kernel.moment(2, 0, ker)
    m3 <- kernel.moment(3, 0, ker)
  }
})

local({
  ## idw
  Z <- idw(longleaf, power=4)
  Z <- idw(longleaf, power=4, se=TRUE)
  ZX <- idw(longleaf, power=4, at="points")
  ZX <- idw(longleaf, power=4, at="points", se=TRUE)
  ## dodgy code blocks in densityVoronoi.R
  A <- adaptive.density(nztrees, nrep=2, f=0.5, counting=TRUE)
  B <- adaptive.density(nztrees, nrep=2, f=0.5, counting=TRUE, fixed=TRUE)
  D <- adaptive.density(nztrees, nrep=2, f=0.5, counting=FALSE)
  E <- adaptive.density(nztrees, nrep=2, f=0.5, counting=FALSE, fixed=TRUE)
  #' adaptive kernel estimation
  d10 <- nndist(nztrees, k=10)
  d10fun <- distfun(nztrees, k=10)
  d10im  <- as.im(d10fun)
  uN <- 2 * runif(npoints(nztrees))
  AA <- densityAdaptiveKernel(nztrees, bw=d10)
  BB <- densityAdaptiveKernel(nztrees, bw=d10, weights=uN)
  DD <- densityAdaptiveKernel(nztrees, bw=d10fun, weights=uN)
  EE <- densityAdaptiveKernel(nztrees, bw=d10im, weights=uN)
})

local({
  ## unnormdensity
  x <- rnorm(20) 
  d0 <- unnormdensity(x, weights=rep(0, 20))
  dneg <- unnormdensity(x, weights=c(-runif(19), 0))
})

local({
  ## cases of 'intensity' etc
  a <- intensity(amacrine, weights=expression(x))
  a <- intensity(split(amacrine), weights=expression(x))
  a <- intensity(split(amacrine), weights=amacrine$x)
  a <- intensity(ppm(amacrine ~ 1))
})
reset.spatstat.options()

#'
#'      tests/diagnostique.R
#'
#'  Diagnostic tools such as diagnose.ppm, qqplot.ppm
#'
#'  $Revision: 1.5 $  $Date: 2019/12/31 03:32:54 $
#'

require(spatstat)
local({
  fit <- ppm(cells ~ x)
  diagE <- diagnose.ppm(fit, type="eem")
  diagI <- diagnose.ppm(fit, type="inverse")
  diagP <- diagnose.ppm(fit, type="Pearson")
  plot(diagE, which="all")
  plot(diagI, which="smooth")
  plot(diagP, which="x")
  plot(diagP, which="marks", plot.neg="discrete")
  plot(diagP, which="marks", plot.neg="contour")
  plot(diagP, which="smooth", srange=c(-5,5))
  plot(diagP, which="smooth", plot.smooth="contour")
  plot(diagP, which="smooth", plot.smooth="image")

  fitS <- ppm(cells ~ x, Strauss(0.08))
  diagES <- diagnose.ppm(fitS, type="eem", clip=FALSE)
  diagIS <- diagnose.ppm(fitS, type="inverse", clip=FALSE)
  diagPS <- diagnose.ppm(fitS, type="Pearson", clip=FALSE)
  plot(diagES, which="marks", plot.neg="imagecontour")
  plot(diagPS, which="marks", plot.neg="discrete")
  plot(diagPS, which="marks", plot.neg="contour")
  plot(diagPS, which="smooth", plot.smooth="image")
  plot(diagPS, which="smooth", plot.smooth="contour")
  plot(diagPS, which="smooth", plot.smooth="persp")
  
  #' infinite reach, not border-corrected
  fut <- ppm(cells ~ x, Softcore(0.5), correction="isotropic")
  diagnose.ppm(fut)

  #' 
  diagPX <- diagnose.ppm(fit, type="Pearson", cumulative=FALSE)
  plot(diagPX, which="y")

  #' simulation based
  e <- envelope(cells, nsim=4, savepatterns=TRUE, savefuns=TRUE)
  Plist <- rpoispp(40, nsim=5)

  qf <- qqplot.ppm(fit, nsim=4, expr=e, plot.it=FALSE)
  print(qf)
  qp <- qqplot.ppm(fit, nsim=5, expr=Plist, fast=FALSE)
  print(qp)
  qp <- qqplot.ppm(fit, nsim=5, expr=expression(rpoispp(40)), plot.it=FALSE)
  print(qp)
  qg <- qqplot.ppm(fit, nsim=5, style="classical", plot.it=FALSE)
  print(qg)
  
  #' lurking.ppm
  #' covariate is numeric vector
  fitx <- ppm(cells ~ x)
  yvals <- coords(as.ppp(quad.ppm(fitx)))[,"y"]
  lurking(fitx, yvals)
  #' covariate is stored but is not used in model
  Z <- as.im(function(x,y){ x+y }, Window(cells))
  fitxx <- ppm(cells ~ x, data=solist(Zed=Z), allcovar=TRUE)
  lurking(fitxx, expression(Zed))
  #' envelope is a ppplist; length < nsim; glmdata=NULL
  fit <- ppm(cells ~ 1)
  stuff <- lurking(fit, expression(x), envelope=Plist, plot.sd=FALSE)
  #' plot.lurk
  plot(stuff, shade=NULL)
})
#'
#'  tests/discarea.R
#'
#'   $Revision: 1.2 $ $Date: 2019/01/20 08:44:50 $
#'

require(spatstat)
local({
  u <- c(0.5,0.5)
  B <- owin(poly=list(x=c(0.3, 0.5, 0.7, 0.4), y=c(0.3, 0.3, 0.6, 0.8)))
  areaGain(u, cells, 0.1, exact=TRUE)
  areaGain(u, cells, 0.1, W=NULL)
  areaGain(u, cells, 0.1, W=B)

  X <- cells[square(0.4)]
  areaLoss(X, 0.1, exact=TRUE)  # -> areaLoss.diri
  areaLoss(X, 0.1, exact=FALSE) # -> areaLoss.grid
  areaLoss.poly(X, 0.1)

  areaLoss(X, 0.1, exact=FALSE, method="distmap")          # -> areaLoss.grid
  areaLoss(X, c(0.1, 0.15), exact=FALSE, method="distmap") # -> areaLoss.grid
})
#'
#'   tests/disconnected.R
#'
#'   disconnected linear networks
#'
#'    $Revision: 1.3 $ $Date: 2018/07/21 03:00:09 $

require(spatstat)

local({

#'   disconnected network
m <- simplenet$m
m[4,5] <- m[5,4] <- m[6,10] <- m[10,6] <- m[4,6] <- m[6,4] <- FALSE
L <- linnet(vertices(simplenet), m)
L
summary(L)
is.connected(L)
Z <- connected(L, what="components")

#' point pattern with no points in one connected component
set.seed(42)
X <- rpoislpp(lambda=function(x,y) { 10 * (x < 0.5)}, L)
B <- lineardirichlet(X)
plot(B)
summary(B)
D <- pairdist(X)
A <- nndist(X)
H <- nnwhich(X)
Y <- rpoislpp(lambda=function(x,y) { 10 * (x < 0.5)}, L)
G <- nncross(X, Y)
J <- crossdist(X, Y)
plot(distfun(X))  # includes evaluation of nncross(what="dist")

#' K functions in disconnected network
K <- linearK(X)
lamX <- intensity(X)
nX <- npoints(X)
KI <- linearKinhom(X, lambda=rep(lamX, nX))
P <- linearpcf(X)
PJ <- linearpcfinhom(X, lambda=rep(lamX, nX))
Y <- X %mark% factor(rep(1:2, nX)[1:nX])
Y1 <- split(Y)[[1]]
Y2 <- split(Y)[[2]]
KY <- linearKcross(Y)
PY <- linearpcfcross(Y)
KYI <- linearKcross.inhom(Y, lambdaI=rep(intensity(Y1), npoints(Y1)),
                       lambdaJ=rep(intensity(Y2), npoints(Y2)))
PYI <- linearpcfcross.inhom(Y, lambdaI=rep(intensity(Y1), npoints(Y1)),
                    lambdaJ=rep(intensity(Y2), npoints(Y2)))

#' internal utilities
K <- ApplyConnected(X, linearK, rule=function(...) list())
})
#'
#'   tests/dominic.R
#'
#'   Additional tests for Dominic Schuhmacher's code
#'
#'   $Revision: 1.3 $  $Date: 2019/10/11 04:33:29 $

require(spatstat)
local({
  X <- runifpoint(10)
  Y <- runifpoint(10)

  d  <- pppdist(X, Y, type="ace", show.rprimal=TRUE)
  a <- matchingdist(d, type="ace")
  b <- matchingdist(d, type="mat")

  d2 <- pppdist(X, Y, type="spa", ccode=FALSE)
  d2 <- pppdist(X, Y, type="spa", ccode=TRUE, auction=FALSE)
  d3 <- pppdist(X, Y, type="mat", ccode=TRUE, auction=FALSE)
  d4 <- pppdist(X[FALSE], Y[FALSE], matching=TRUE, type="spa")
  d4 <- pppdist(X[FALSE], Y[FALSE], matching=FALSE, type="spa")
  d4 <- pppdist(X[FALSE], Y[FALSE], matching=TRUE, type="ace")
  d4 <- pppdist(X[FALSE], Y[FALSE], matching=FALSE, type="ace")

  m  <- pppdist.mat(X, Y, q=Inf, cutoff=0.001)
  m2 <- pppdist.mat(X[FALSE], Y[FALSE], q=Inf, cutoff=0.001)
  m3 <- pppdist.mat(X[FALSE], Y[FALSE], q=2, cutoff=0.001)

})



#'
#'    tests/deepeepee.R
#'
#'    Tests for determinantal point process models
#' 
#'    $Revision: 1.6 $ $Date: 2020/01/10 03:10:21 $

require(spatstat)
local({
  #' simulate.dppm
  jpines <- residualspaper$Fig1
  fit <- dppm(jpines ~ 1, dppGauss)
  set.seed(10981)
  simulate(fit, W=square(5))
  #' simulate.detpointprocfamily - code blocks
  model <- dppGauss(lambda=100, alpha=.05, d=2)
  simulate(model, seed=1999, correction="border")
  u <- is.stationary(model)
  #' other methods for dppm
  kay <- Kmodel(fit)
  gee <- pcfmodel(fit)
  lam <- intensity(fit)
  arr <- reach(fit)
  pah <- parameters(fit)
  
  #' dppeigen code blocks
  mod <- dppMatern(lambda=2, alpha=0.01, nu=1, d=2)
  uT <- dppeigen(mod, trunc=1.1,  Wscale=c(1,1), stationary=TRUE)
  uF <- dppeigen(mod, trunc=1.1,  Wscale=c(1,1), stationary=FALSE)
  vT <- dppeigen(mod, trunc=0.98, Wscale=c(1,1), stationary=TRUE)
  vF <- dppeigen(mod, trunc=0.98, Wscale=c(1,1), stationary=FALSE)
})
#'
#'   tests/duplicity.R
#'
#'  Tests of duplicated/multiplicity code
#'
#' $Revision: 1.7 $ $Date: 2019/12/06 02:41:32 $

require(spatstat)
local({
   X <- ppp(c(1,1,0.5,1), c(2,2,1,2), window=square(3), check=FALSE)
   Y <- X %mark% factor(letters[c(3,2,4,3)])
   ZC <- X %mark% letters[c(3,2,4,3)]
   ZM <- Y %mark% matrix(c(3,2,4,3), 4, 2)
   ZD <- Y %mark% as.data.frame(marks(ZM))

   #' multiplicity
   m <- multiplicity(X)
   mf <- multiplicity(Y)
   mm <- multiplicity(ZM)
   mz <- multiplicity(ZD)
   mc <- multiplicity(ZC)
   ## default method
   kk <- c(1,2,3,1,1,2)
   mk <- multiplicity(kk)
   ml <- multiplicity(list(sin, cos, tan)[kk])
   mc <- multiplicity(c("sin", "cos", "tan")[kk])
   if(!identical(ml, mk))
     stop("multiplicity.default(<list>) disagrees with multiplicityNumeric")
   if(!identical(mc, mk))
     stop("multiplicity(<character>) disagrees with multiplicity(<numeric>)")
   ## data frame method
   df <- data.frame(x=c(1:4, 1,3,2,4, 0,0, 3,4),
                    y=factor(rep(letters[1:4], 3)))
   md <- multiplicity(df)

   ## uniquemap.ppp
   checkum <- function(X, blurb) {
     a <- uniquemap(X)
     if(any(a > seq_along(a)))
       stop(paste("uniquemap", blurb,
                  "does not respect sequential ordering"))
     return(invisible(NULL))
   }
   checkum(X, "<unmarked point pattern>")
   checkum(Y, "<multitype point pattern>")
   checkum(ZC, "<point pattern with character marks>")
   checkum(ZM, "<point pattern with matrix of marks>")
   checkum(ZD, "<point pattern with several columns of marks>")

   ## uniquemap.data.frame
   dfbase <- as.data.frame(replicate(3, sample(1:20, 10), simplify=FALSE))
   df <- dfbase[sample(1:10, 30, replace=TRUE), , drop=FALSE]
   #' faster algorithm for numeric values
   checkum(df, "<numeric data frame>")
   a <- uniquemap(df)
   #' general algorithm using 'duplicated' and 'match'
   dfletters <- as.data.frame(matrix(letters[as.matrix(df)], nrow=nrow(df)))
   checkum(dfletters, "<character data frame>")
   b <- uniquemap(dfletters)
   if(!isTRUE(all.equal(a,b)))
     stop("inconsistency between algorithms in uniquemap.data.frame")

   ## uniquemap.matrix
   M0 <- matrix(1:12, 3, 4)
   ii <- sample(1:3, 5, replace=TRUE)
   M4 <- M0[ii, , drop=FALSE]
   checkum(M4, "<integer matrix>")
   u4 <- uniquemap(M4)
   C4 <- matrix(letters[M4], 5, 4)
   uc4 <- uniquemap(C4)
   checkum(C4, "<character matrix>")
   if(!isTRUE(all.equal(u4, uc4)))
     stop("Inconsistency between algorithms in uniquemap.matrix")
   
   ## uniquemap.default
   a <- letters[c(1, 1:4, 3:2)]
   checkum(a, "<character>")
   checkum(as.list(a), "<list>")
   u1 <- uniquemap(a)
   u2 <- uniquemap(as.list(a))
   if(!isTRUE(all.equal(u1, u2)))
     stop("Inconsistency between algorithms in uniquemap.default")
})
