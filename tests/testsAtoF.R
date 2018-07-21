#'  tests/aucroc.R
#'
#'  AUC and ROC code
#'
#'  $Revision: 1.1 $ $Date: 2018/05/13 23:56:20 $

require(spatstat)
local({
  fit <- kppm(redwood ~ I(y-x))
  a <- roc(fit)
  b <- auc(fit)
  d <- roc(spiders, "x")
  e <- auc(spiders, "y")
  fut <- lppm(spiders ~ I(y-x))
  f <- roc(fut)
  g <- roc(fut)
})
## badwindowcheck.R
## $Revision: 1.2 $  $Date: 2014/01/27 07:18:41 $
##

require(spatstat)
local({
  ## Simple example of self-crossing polygon
  x <- read.table("selfcross.txt", header=TRUE)
  ## Auto-repair
  w <- owin(poly=x)

  ## Real data involving various quirks
  b <- read.table("badwindow.txt", header=TRUE)
  b <- split(b, factor(b$i))
  b <- lapply(b, function(z) { as.list(z[,-3]) })
  ## make owin without checking
  W <- owin(poly=b, check=FALSE)
  ## Apply stringent checks
  owinpolycheck(W,verbose=FALSE)
  ## Auto-repair
  W2 <- owin(poly=b)
})




## tests/cdf.test.R
require(spatstat)
local({
  ## (1) check cdf.test with strange data
  ## Marked point patterns with some marks not represented
  AC <- split(ants, un=FALSE)$Cataglyphis
  AM <- split(ants, un=FALSE)$Messor
  DM <- distmap(AM)
  ## should produce a warning, rather than a crash:
  cdf.test(AC, DM)
  ## should be OK:
  cdf.test(unmark(AC), DM)
  cdf.test(unmark(AC), DM, "cvm")
  cdf.test(unmark(AC), DM, "ad")
  ## other code blocks
  cdf.test(finpines, "x")

  ## (2) linear networks
  set.seed(42)
  X <- runiflpp(20, simplenet)
  cdf.test(X, "x")
  cdf.test(X, "x", "cvm")
  fit <- lppm(X ~1)
  cdf.test(fit, "y")
  cdf.test(fit, "y", "cvm")
  cdf.test(fit, "y", "ad")
  ## marked
  cdf.test(chicago, "y")

  ## (3) Monte Carlo test for Gibbs model
  fit <- ppm(cells ~ 1, Strauss(0.07))
  cdf.test(fit, "x", nsim=9)

  ## cdf.test.slrm
  fut <- slrm(japanesepines ~ x + y)
  Z <- distmap(japanesepines)
  cdf.test(fut, Z)
})



##  tests/closeshave.R
## check 'closepairs/crosspairs' code
## validity and memory allocation
## $Revision: 1.13 $ $Date: 2018/06/07 05:55:00 $

local({
  r <- 0.12
  close.all <- closepairs(redwood, r)
  close.ij <- closepairs(redwood, r, what="indices")
  close.ijd <- closepairs(redwood, r, what="ijd")
  close.every <- closepairs(redwood, r, what="all", distinct=FALSE)
  stopifnot(identical(close.ij, close.all[c("i","j")]))
  stopifnot(identical(close.ijd, close.all[c("i","j","d")]))

  Y <- split(amacrine)
  on <- Y$on
  off <- Y$off
  cross.all <- crosspairs(on, off, r)
  cross.ij <- crosspairs(on, off, r, what="indices")
  cross.ijd <- crosspairs(on, off, r, what="ijd")
  cross.every <- crosspairs(on, off, r, what="all", distinct=FALSE)
  stopifnot(identical(cross.ij, cross.all[c("i","j")]))
  stopifnot(identical(cross.ijd, cross.all[c("i","j","d")]))

  # closethresh vs closepairs: EXACT agreement
  thresh <- 0.08
  clt <- closethresh(redwood, r, thresh)
  cl <- with(closepairs(redwood, r),
             list(i=i, j=j, th = (d <= thresh)))
  if(!identical(cl, clt))
    stop("closepairs and closethresh disagree")

  reordered <- function(a) {
    o <- with(a, order(i,j))
    as.list(as.data.frame(a)[o,,drop=FALSE])
  }
  samesame <- function(a, b) {
    identical(reordered(a), reordered(b))
  }
  
  ## ...............................................
  #' compare with older, slower code
  op <- spatstat.options(closepairs.newcode=FALSE,
                         crosspairs.newcode=FALSE)
  ## ...............................................
  old.close.ij <- closepairs(redwood, r, what="indices")
  old.cross.ij <- crosspairs(on, off, r, what="indices")
  stopifnot(samesame(close.ij, old.close.ij))
  stopifnot(samesame(cross.ij, old.cross.ij))
  # execute only:
  old.close.every <- closepairs(redwood, r, what="all", distinct=FALSE)
  old.cross.every <- crosspairs(on, off, r, what="all", distinct=FALSE)
  ## ...............................................
  spatstat.options(op)
  ## ...............................................

  ## ...............................................
  #' alternative code - execution only
  op <- spatstat.options(closepairs.newcode=FALSE,
                         closepairs.altcode=TRUE)
  alt.close.ij <- closepairs(redwood, r, what="indices")
  alt.close.ijd <- closepairs(redwood, r, what="ijd")
  alt.close.all <- closepairs(redwood, r, what="all")
  spatstat.options(op)
  ## ...............................................
  
  # Rasmus' example
  R <- 0.04
  U <- as.ppp(gridcenters(owin(), 50, 50), W=owin())
  cp <- crosspairs(U, U, R)
  G <- matrix(0, npoints(U), npoints(U))
  G[cbind(cp$i, cp$j)] <- 1
  if(!isSymmetric(G))
    stop("crosspairs is not symmetric in Rasmus example")

  #' periodic distance
  pclose <- function(X, R, method=c("raw", "C")) {
    method <- match.arg(method)
    switch(method,
           raw = {
             D <- pairdist(X, periodic=TRUE)
             diag(D) <- Inf
             result <- which(D <= R, arr.ind=TRUE)
           },
           C = {
             result <- closepairs(X, R, periodic=TRUE, what="indices")
           })
    result <- as.data.frame(result)
    colnames(result) <- c("i","j")
    return(result)
  }
  #' pick a threshold value which avoids GCC bug 323
  RR <- 0.193
  A <- pclose(redwood, RR, "raw")
  B <- pclose(redwood, RR, "C")
  if(!samesame(A,B))
    stop("closepairs.ppp(periodic=TRUE) gives wrong answer")

  #' other functions that don't have a help file
  niets <- crosspairquad(quadscheme(cells), 0.1)

  #' other code blocks
  u <- closepairs(cells, 0.09, periodic=TRUE, what="all")
  v <- closepairs(cells, 0.07, twice=FALSE, neat=TRUE)
  #' tight cluster - guess count does not work
  Xc <- runifpoint(100, square(0.01))
  Window(Xc) <- square(1)
  z <- closepairs(Xc, 0.02, what="indices", distinct=FALSE)
  z <- closepairs(Xc, 0.02, what="ijd",     distinct=FALSE)
  z <- closepairs(Xc, 0.02, what="all",     distinct=FALSE)
  #' same task, older code
  aop <- spatstat.options(closepairs.newcode=FALSE)
  z <- closepairs(Xc, 0.02, what="indices", distinct=FALSE)
  z <- closepairs(Xc, 0.02, what="ijd",     distinct=FALSE)
  z <- closepairs(Xc, 0.02, what="all",     distinct=FALSE)
  spatstat.options(aop)
})

reset.spatstat.options()
## tests/colour.R
##
##  Colour value manipulation and colour maps
##
## $Revision: 1.3 $ $Date: 2018/06/08 13:33:39 $
##

require(spatstat)

local({
   f <- function(n) grey(seq(0,1,length=n))
   z <- to.grey(f)

   a <- colourmap(rainbow(12), range=as.Date(c("2018-01-01", "2018-12-31")))
   print(a)
   print(summary(a))

   b <- colourmap(rainbow(12), inputs=month.name)
   plot(b, vertical=FALSE)
   plot(b, vertical=TRUE)
   
})
#'
#'   tests/contrib.R
#'
#'   Tests for user-contributed code in spatstat
#'
#'   $Revision: 1.1 $  $Date: 2018/07/01 04:48:25 $

require(spatstat)
local({
  #' Jinhom
  #' Marie-Colette van Lieshout and Ottmar Cronie
  X <- redwood3
  fit <- ppm(X ~ polynom(x,y,2))
  lam <- predict(fit)
  lamX <- fitted(fit, dataonly=TRUE)
  lmin <- 0.9 * min(lam)
  g1 <- Ginhom(X, lambda=fit, update=TRUE)
  g2 <- Ginhom(X, lambda=fit, update=FALSE, lmin = lmin)
  g3 <- Ginhom(X, lambda=lam,  lmin=lmin)
  g4 <- Ginhom(X, lambda=lamX, lmin=lmin)
  f1 <- Finhom(X, lambda=fit, update=TRUE)
  f2 <- Finhom(X, lambda=fit, update=FALSE)
  f3 <- Finhom(X, lambda=lam,  lmin=lmin)
})
# tests/correctC.R
# check for agreement between C and interpreted code
# for interpoint distances etc.
# $Revision: 1.4 $ $Date: 2015/12/29 08:54:49 $

require(spatstat)

local({
  eps <- .Machine$double.eps * 4

  # pairdist.ppp
  X <- rpoispp(42)
  dC <- pairdist(X, method="C")
  dR <- pairdist(X, method="interpreted")
  if(any(abs(dC - dR) > eps))
    stop("Algorithms for pairdist() do not agree")

  dC <- pairdist(X, periodic=TRUE, method="C")
  dR <- pairdist(X, periodic=TRUE, method="interpreted")
  if(any(abs(dC - dR) > eps))
    stop("Algorithms for pairdist(periodic=TRUE) do not agree")

  # crossdist.ppp
  Y <- rpoispp(42)
  dC <- crossdist(X, Y, method="C")
  dR <- crossdist(X, Y, method="interpreted")
  if(any(abs(dC - dR) > eps))
    stop("Algorithms for crossdist() do not agree")

  dC <- crossdist(X, Y, periodic=TRUE, method="C")
  dR <- crossdist(X, Y, periodic=TRUE, method="interpreted")
  if(any(abs(dC - dR) > eps))
    stop("Algorithms for crossdist(periodic=TRUE) do not agree")

  # nndist.ppp
  nnC <- nndist(X, method="C")
  nnI <- nndist(X, method="interpreted")
  if(any(abs(nnC - nnI) > eps))
    stop("Algorithms for nndist() do not agree")

  nn3C <- nndist(X, k=3, method="C")
  nn3I <- nndist(X, k=3, method="interpreted")
  if(any(abs(nn3C - nn3I) > eps))
    stop("Algorithms for nndist(k=3) do not agree")

  # nnwhich.ppp
  nwC <- nnwhich(X, method="C")
  nwI <- nnwhich(X, method="interpreted")
  if(any(nwC != nwI))
    stop("Algorithms for nnwhich() do not agree")

  nw3C <- nnwhich(X, k=3, method="C")
  nw3I <- nnwhich(X, k=3, method="interpreted")
  if(any(nw3C != nw3I))
    stop("Algorithms for nnwhich(k=3) do not agree")

  # whist
  set.seed(98123)
  x <- runif(1000)
  w <- sample(1:5, 1000, replace=TRUE)
  b <- seq(0,1,length=101)
  op <- spatstat.options(Cwhist=TRUE)
  aT <- whist(x,b,w)
  spatstat.options(Cwhist=FALSE)
  aF <- whist(x,b,w)
  if(!all(aT == aF))
    stop("Algorithms for whist disagree")
  spatstat.options(op)
})

reset.spatstat.options()

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
#'
#'  $Revision: 1.22 $  $Date: 2018/07/12 02:14:08 $
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

  V <- diag(c(0.05^2, 0.07^2))
  tryit(varcov=V)
  tryit(varcov=V, diggle=TRUE)
  tryit(varcov=V, weights=expression(x))
  tryit(varcov=V, weights=expression(x), diggle=TRUE)
  Z <- distmap(runifpoint(5, Window(cells)))
  tryit(0.05, weights=Z)
  tryit(0.05, weights=Z, diggle=TRUE)

  wdf <- data.frame(a=1:42,b=42:1)
  tryit(0.05, weights=wdf, do.fun=FALSE)
  tryit(varcov=V, weights=wdf, do.fun=FALSE)
  tryit(varcov=V, weights=expression(cbind(x,y)), do.fun=FALSE)

  crossit <- function(..., sigma=NULL) {
    U <- runifpoint(20, Window(cells))
    a <- densitycrossEngine(cells, U, sigma=sigma, ...)
    a <- densitycrossEngine(cells, U, sigma=sigma, ..., diggle=TRUE)
    invisible(NULL)
  }
  crossit(varcov=V, weights=cells$x)
  crossit(varcov=V, weights=wdf)

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

  lam <- density(redwood)
  K <- Kinhom(redwood, lam)
  
  lamX <- density(redwood, at="points")
  KX <- Kinhom(redwood, lamX)

  ## test all code cases of new 'relrisk.ppp' algorithm
  pants <- function(..., X=ants) {
    a <- relrisk(X, sigma=100, se=TRUE, ...)
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

  ## more than 2 types
  pants(X=sporophores)
  pants(X=sporophores, at="points")
  pants(X=sporophores, relative=TRUE, at="points")

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
  stroke <- function(...) {
    Z <- Smooth(longleaf, ..., at="pixels")
    Z <- Smooth(longleaf, ..., at="points")
    f <- Smoothfun(longleaf, ...)
    f(120, 80)
    return(invisible(NULL))
  }
  stroke()
  stroke(5, diggle=TRUE)
  stroke(1e-6) # generates warning about small bandwidth
  stroke(varcov=diag(c(25, 36)))
  stroke(5, weights=runif(npoints(longleaf)))
  stroke(5, weights=expression(x))
  strike <- function(...) {
    Z <- Smooth(finpines, ..., at="pixels")
    Z <- Smooth(finpines, ..., at="points")
    f <- Smoothfun(finpines, ...)
    f(4, 1)
    return(invisible(NULL))
  }
  strike()
  strike(varcov=diag(c(1.2, 2.1)))
  strike(1.5, weights=runif(npoints(finpines)))
  strike(1.5, weights=expression(y))

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
  A <- second.moment.calc(redwood, 0.1, what="all", debug=TRUE)
  B <- second.moment.calc(redwood, varcov=diag(c(0.1,0.1)^2), what="all")
  PR <- pixellate(redwood)
  DR <- second.moment.calc(list(PR, PR), 0.1, debug=TRUE,
                             npts=npoints(redwood), obswin=Window(redwood))
})

local({
  #' bandwidth selection
  op <- spatstat.options(n.bandwidth=8)
#  bw.relrisk(urkiola, hmax=20) is tested in man/bw.relrisk.Rd
  bw.relrisk(urkiola, hmax=20, method="leastsquares")
  bw.relrisk(urkiola, hmax=20, method="weightedleastsquares")
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

reset.spatstat.options()

#'
#'      tests/diagnostique.R
#'
#'  Diagnostic tools such as diagnose.ppm, qqplot.ppm
#'
#'  $Revision: 1.1 $  $Date: 2018/05/22 11:57:12 $
#'

require(spatstat)
local({
  fit <- ppm(cells ~ x)
  e <- envelope(cells, nsim=4, savepatterns=TRUE, savefuns=TRUE)
  qf <- qqplot.ppm(fit, nsim=4, expr=e, plot.it=FALSE)
  print(qf)
  qg <- qqplot.ppm(fit, nsim=5, style="classical", plot.it=FALSE)
})
#'
#'  tests/discarea.R
#'
#'   $Revision: 1.1 $ $Date: 2016/03/28 09:16:03 $
#'

require(spatstat)
local({
  u <- c(0.5,0.5)
  B <- owin(poly=list(x=c(0.3, 0.5, 0.7, 0.4), y=c(0.3, 0.3, 0.6, 0.8)))
  areaGain(u, cells, 0.1, exact=TRUE)
  areaGain(u, cells, 0.1, W=NULL)
  areaGain(u, cells, 0.1, W=B)

  areaLoss(cells[square(0.4)], 0.1, exact=TRUE)
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
#'   $Revision: 1.1 $  $Date: 2018/04/29 11:01:30 $

require(spatstat)
local({
  X <- runifpoint(10)
  Y <- runifpoint(10)
  d <- pppdist(X, Y, type="ace", show.rprimal=TRUE)
  m <- pppdist.mat(X, Y, q=Inf, cutoff=0.001)
})
#  tests/emptymarks.R
#
# test cases where there are no (rows or columns of) marks
#
#  $Revision: 1.3 $ $Date: 2015/12/29 08:54:49 $

require(spatstat)
local({
  n <- npoints(cells)
  df <- data.frame(x=1:n, y=factor(sample(letters, n, replace=TRUE)))
  nocolumns <- c(FALSE, FALSE)
  norows <- rep(FALSE, n)
  X <- cells
  marks(X) <- df
  marks(X) <- df[,1]
  marks(X) <- df[,nocolumns]
  Z <- Y <- X[integer(0)]
  marks(Y) <- df[norows,]
  stopifnot(is.marked(Y))
  marks(Z) <- df[norows,nocolumns]
  stopifnot(!is.marked(Z))
})
#
#  tests/envelopes.R
#
#  Test validity of envelope data
#
#  $Revision: 1.11 $  $Date: 2018/07/12 02:14:13 $
#

require(spatstat)

local({
checktheo <- function(fit) {
  fitname <- deparse(substitute(fit))
  en <- envelope(fit, nsim=4, verbose=FALSE, nrep=1e3)
  nama <- names(en)
  expecttheo <- is.poisson(fit) && is.stationary(fit)
  context <- paste("Envelope of", fitname)
  if(expecttheo) {
    if(!("theo" %in% nama))
      stop(paste(context, "did not contain", sQuote("theo")))
    if("mmean" %in% nama)
      stop(paste(context, "unexpectedly contained", sQuote("mmean")))
  } else {
    if("theo" %in% nama)
      stop(paste(context, "unexpectedly contained", sQuote("theo")))
    if(!("mmean" %in% nama))
      stop(paste(context, "did not contain", sQuote("mmean")))
  }
  cat(paste(context, "has correct format\n"))
}
  
checktheo(ppm(cells))
checktheo(ppm(cells ~x))
checktheo(ppm(cells ~1, Strauss(0.1)))

# check envelope calls from 'alltypes'
a <- alltypes(demopat, Kcross, nsim=4, envelope=TRUE)
b <- alltypes(demopat, Kcross, nsim=4, envelope=TRUE, global=TRUE)

# check 'transform' idioms
A <- envelope(cells, Kest, nsim=4, transform=expression(. - .x))
B <- envelope(cells, Kest, nsim=4, transform=expression(sqrt(./pi) - .x))

#' check savefuns/savepatterns with global 
fit <- ppm(cells~x)
Ef <- envelope(fit, Kest, nsim=4, savefuns=TRUE, global=TRUE)
Ep <- envelope(fit, Kest, nsim=4, savepatterns=TRUE, global=TRUE)

#' check handling of 'dangerous' cases
fut <- ppm(redwood ~ x)
Ek <- envelope(fut, Kinhom, update=FALSE, nsim=4)

# check conditional simulation
e1 <- envelope(cells, Kest, nsim=4, fix.n=TRUE)
e2 <- envelope(amacrine, Kest, nsim=4, fix.n=TRUE)
e3 <- envelope(amacrine, Kcross, nsim=4, fix.marks=TRUE)
e4 <- envelope(finpines, Kest, nsim=4, fix.n=TRUE) # multiple columns of marks
e5 <- envelope(finpines, Kest, nsim=4, fix.marks=TRUE)
fit <- ppm(japanesepines ~ 1, Strauss(0.04))
e6 <- envelope(fit, Kest, nsim=4, fix.n=TRUE)
fit2 <- ppm(amacrine ~ 1, Strauss(0.03))
e7 <- envelope(fit2, Gcross, nsim=4, fix.marks=TRUE)

# check pooling of envelopes in global case
E1 <- envelope(cells, Kest, nsim=5, savefuns=TRUE, global=TRUE)
E2 <- envelope(cells, Kest, nsim=12, savefuns=TRUE, global=TRUE)
p12 <- pool(E1, E2)
p12 <- pool(E1, E2, savefuns=TRUE)
F1 <- envelope(cells, Kest, nsim=5,
               savefuns=TRUE, savepatterns=TRUE, global=TRUE)
F2 <- envelope(cells, Kest, nsim=12,
               savefuns=TRUE, savepatterns=TRUE, global=TRUE)
p12 <- pool(F1, F2)
p12 <- pool(F1, F2, savefuns=TRUE, savepatterns=TRUE)
E1r <- envelope(cells, Kest, nsim=5, savefuns=TRUE, global=TRUE,
                ginterval=c(0.05, 0.15))
E2r <- envelope(cells, Kest, nsim=12, savefuns=TRUE, global=TRUE,
                ginterval=c(0.05, 0.15))
p12r <- pool(E1r, E2r)
})

local({
  #' as.data.frame.envelope
  Nsim <- 5
  E <- envelope(cells, nsim=Nsim, savefuns=TRUE)
  A <- as.data.frame(E)
  B <- as.data.frame(E, simfuns=TRUE)
  stopifnot(ncol(B) - ncol(A) == Nsim)
})

local({
  #' cases not covered elsewhere
  A <- envelope(cells, nsim=5, alternative="less",
                do.pwrong=TRUE, use.theory=FALSE,
                savepatterns=TRUE, savefuns=TRUE)
  print(A)
  B <- envelope(A, nsim=5, savefuns=TRUE)
  D <- envelope(cells, "Lest", nsim=5)
  fit <- ppm(cells ~ 1, Strauss(0.07))
  U <- envelope(fit, nsim=3, simulate=expression(runifpoint(20)))
  #'  envelopes based on sample variance
  E <- envelope(cells, nsim=8, VARIANCE=TRUE)
  G <- envelope(cells, nsim=8, VARIANCE=TRUE,
                use.theory=FALSE, do.pwrong=TRUE)
  print(G)
  #' summary method
  summary(E)
  summary(envelope(cells, nsim=5, simulate=expression(runifpoint(42))))
  #' weights argument
  H1 <- envelope(cells, nsim=4, weights=npoints, savefuns=TRUE)
  H2 <- envelope(cells, nsim=4, weights=npoints, savefuns=TRUE)
  J <- envelope(cells, nsim=4, weights=npoints, VARIANCE=TRUE)
  #' pooling with weights
  H <- pool(H1, H2)
  #' pooling envelopes with non-identical attributes
  H0 <- envelope(cells, nsim=4, savefuns=TRUE)
  HH <- pool(H0, H1)
  #' undocumented/secret
  K <- envelope(cells, nsim=4, saveresultof=npoints, collectrubbish=TRUE)
  #' so secret I've even forgotten how to do it
  M <- envelope(cells, nsim=4, internal=list(eject="patterns"))
})

local({
  #' envelope computations in other functions
  P <- lurking(cells, expression(x), envelope=TRUE, nsim=9)
  print(P)
  #' re-using envelope objects in other functions
  A <- envelope(cells, nsim=9, savepatterns=TRUE, savefuns=TRUE)
  S <- lurking(cells, expression(x), envelope=A, nsim=9)
  #' envelope.envelope
  B <- envelope(cells, nsim=5, savepatterns=TRUE, savefuns=FALSE)
  envelope(B)
})
#
#    tests/factorbugs.R
#
# check for various bugs related to factor conversions
#
#    $Revision: 1.3 $  $Date: 2015/12/29 08:54:49 $
#
require(spatstat)
local({
  # make a factor image
  m <- factor(rep(letters[1:4], 4))
  Z <- im(m, xcol=1:4, yrow=1:4)
  # make a point pattern
  set.seed(42)
  X <- runifpoint(20, win=as.owin(Z))
  # look up the image at the points of X
  # (a) internal
  ans1 <- lookup.im(Z, X$x, X$y)
  stopifnot(is.factor(ans1))
  # (b) user level
  ans2 <- Z[X]
  stopifnot(is.factor(ans2))
  # (c) turn the image into a tessellation
  #  and apply quadratcount
  V <- tess(image = Z)
  quadratcount(X, tess=V)
  # (d) pad image
  Y <- padimage(Z, factor("b", levels=levels(Z)))
  stopifnot(Y$type == "factor")
  U <- padimage(Z, "b")
  stopifnot(U$type == "factor")
})


#
#    tests/fastgeyer.R
#
# checks validity of fast C implementation of Geyer interaction
#
#    $Revision: 1.3 $  $Date: 2015/12/29 08:54:49 $
#
require(spatstat)
local({
  X <- redwood
  Q <- quadscheme(X)
  U <- union.quad(Q)
  EP <- equalpairs.quad(Q)
  G <- Geyer(0.11, 2)
# The value r=0.11 is chosen to avoid hardware numerical effects (gcc bug 323).
# It avoids being close any value of pairdist(redwood).
# The nearest such values are 0.1077.. and 0.1131..
# By contrast if r = 0.1 there are values differing from 0.1 by 3e-17
  a <- pairsat.family$eval(X,U,EP,G$pot,G$par,"border")
  b <-          G$fasteval(X,U,EP,G$pot,G$par,"border")
  if(!all(a==b))
    stop("Results of Geyer()$fasteval and pairsat.family$eval do not match")
# ...
# and again for a non-integer value of 'sat'
# (spotted by Thordis Linda Thorarinsdottir)  
  G <- Geyer(0.11, 2.5)
  a <- pairsat.family$eval(X,U,EP,G$pot,G$par,"border")
  b <-          G$fasteval(X,U,EP,G$pot,G$par,"border")
  if(!all(a==b))
    stop("Results of Geyer()$fasteval and pairsat.family$eval do not match when sat is not an integer")
# and again for sat < 1
# (spotted by Rolf)  
  G <- Geyer(0.11, 0.5)
  a <- pairsat.family$eval(X,U,EP,G$pot,G$par,"border")
  b <-          G$fasteval(X,U,EP,G$pot,G$par,"border")
  if(!all(a==b))
    stop("Results of Geyer()$fasteval and pairsat.family$eval do not match when sat < 1")
})

#
#  tests/fastK.R
#
# check fast and slow code for Kest
#       and options not tested elsewhere
#
#   $Revision: 1.4 $   $Date: 2018/04/13 04:01:25 $
#
require(spatstat)
local({
  ## fast code
  Kb <- Kest(cells, nlarge=0)
  Ku <- Kest(cells, correction="none")
  Kbu <- Kest(cells, correction=c("none", "border"))
  ## slow code, full set of corrections, sqrt transformation, ratios
  Ldd <- Lest(unmark(demopat), correction="all", var.approx=TRUE, ratio=TRUE)
  ## Lotwick-Silverman var approx (rectangular window)
  Loo <- Lest(cells, correction="all", var.approx=TRUE, ratio=TRUE)
  ## Code for large dataset
  nbig <- .Machine$integer.max
  if(!is.null(nbig)) {
    nn <- ceiling(sqrt(nbig))
    if(nn < 1e6) Kbig <- Kest(runifpoint(nn),
                              correction=c("border", "bord.modif", "none"),
                              ratio=TRUE)
  }
  
  ## Kinhom
  lam <- density(cells, at="points", leaveoneout=TRUE)
  ## fast code
  Kib <- Kinhom(cells, lam, nlarge=0)
  Kiu <- Kest(cells, lam, correction="none")
  Kibu <- Kest(cells, lam, correction=c("none", "border"))
  ## slow code
  Lidd <- Linhom(unmark(demopat), sigma=bw.scott)
})


#' tests/formuli.R
#'
#'  Test machinery for manipulating formulae
#' 
#' $Revision: 1.4 $  $Date: 2017/02/20 07:35:47 $

require(spatstat)
local({

  ff <- function(A, deletevar, B) {
    D <- reduceformula(A, deletevar)
    if(!spatstat.utils::identical.formulae(D, B)) {
      AD <- as.expression(substitute(reduceformula(A,d),
                                     list(A=A, d=deletevar)))
      stop(paste(AD, "\n\tyields ", spatstat.utils::pasteFormula(D),
                 " instead of ", spatstat.utils::pasteFormula(B)),
           call.=FALSE)
    }
    invisible(NULL)
  }

  ff(~ x + z, "x", ~z)

  ff(y ~ x + z, "x", y~z)

  ff(~ I(x^2) + z, "x",  ~z)

  ff(y ~ poly(x,2) + poly(z,3), "x", y ~poly(z,3))

})



#
#  tests/func.R
#
#   $Revision: 1.3 $   $Date: 2016/06/10 15:04:08 $
#
#  Tests of 'funxy' infrastructure etc

require(spatstat)
local({
  W <- square(1)
  f1a <- function(x, y) sqrt(x^2 + y^2)
  f1b <- function(x, y) { sqrt(x^2 + y^2) }
  f2a <- function(x, y) sin(x)
  f2b <- function(x, y) { sin(x) } 
  f3a <- function(x, y) sin(x) + cos(x) 
  f3b <- function(x, y) { sin(x) + cos(x) } 
  f4a <- function(x, y) { z <- x + y ; z }
  f4b <- function(x, y) { x + y } 
  F1a <- funxy(f1a, W)
  F1b <- funxy(f1b, W)
  F2a <- funxy(f2a, W)
  F2b <- funxy(f2b, W)
  F3a <- funxy(f3a, W)
  F3b <- funxy(f3b, W)
  F4a <- funxy(f4a, W)
  F4b <- funxy(f4b, W)
  stopifnot(identical(F1a(cells), F1b(cells)))
  stopifnot(identical(F2a(cells), F2b(cells)))
  stopifnot(identical(F3a(cells), F3b(cells)))
  stopifnot(identical(F4a(cells), F4b(cells)))
})




##  
##     tests/funnymarks.R
##
## tests involving strange mark values
## $Revision: 1.5 $ $Date: 2018/06/14 01:50:19 $

require(spatstat)
local({
  ## ppm() where mark levels contain illegal characters
  hyphenated <- c("a", "not-a")
  spaced <- c("U", "non U")
  suffixed <- c("a+", "a*")
  charred <- c("+", "*")

  irad <- matrix(0.1, 2,2)
  hrad <- matrix(0.005, 2, 2)

  tryit <- function(types, X, irad, hrad) { 
    levels(marks(X)) <- types
    fit <- ppm(X ~marks + polynom(x,y,2),
               MultiStraussHard(types=types,iradii=irad,hradii=hrad))
    print(fit)
    print(coef(fit))
    val <- fitted(fit)
    pred <- predict(fit)
    return(invisible(NULL))
  }

  tryit(hyphenated, amacrine, irad, hrad)
  tryit(spaced, amacrine, irad, hrad)
  tryit(suffixed, amacrine, irad, hrad)
  tryit(charred, amacrine, irad, hrad)

  ## marks which are dates
  X <- cells
  n <- npoints(X)
  endoftime <- rep(ISOdate(2001,1,1), n)
  eotDate   <- rep(as.Date("2001-01-01"), n)
  markformat(endoftime)
  markformat(eotDate)
  marks(X) <- endoftime
  print(X)
  Y <- X %mark% data.frame(id=1:42, date=endoftime, dd=eotDate)
  print(Y)
  md <- markformat(endoftime)
  
  ## mark formats
  Z <- Y
  marks(Z) <- marks(Z)[1,,drop=FALSE]
  ms <- markformat(solist(cells, redwood))
  marks(Z) <- factor(1:npoints(Z))
  marks(Z)[12] <- NA
  mz <- is.multitype(Z)
  cZ <- coerce.marks.numeric(Z)
  stopifnot(is.multitype(cells %mark% data.frame(a=factor(1:npoints(cells)))))

  a <- numeric.columns(finpines)
  b1 <- numeric.columns(amacrine)
  b2 <- coerce.marks.numeric(amacrine)
  d <- numeric.columns(cells)
  f <- numeric.columns(longleaf)
})
# 
#    tests/fvproblems.R
#
#    $Revision: 1.7 $  $Date: 2016/03/08 00:26:23 $
#

require(spatstat)

# This appears in the workshop notes
# Problem detected by Martin Bratschi

local({
  Jdif <- function(X, ..., i) {
    Jidot <- Jdot(X, ..., i=i)
    J <- Jest(X, ...)
    dif <- eval.fv(Jidot - J)
    return(dif)
  }
  Z <- Jdif(amacrine, i="on")
})
#
#  Test mathlegend code
#
local({
  K <- Kest(cells)
  plot(K)
  plot(K, . ~ r)
  plot(K, . - theo ~ r)
  plot(K, sqrt(./pi)  ~ r)
  plot(K, cbind(iso, theo) ~ r)
  plot(K, cbind(iso, theo) - theo ~ r)
  plot(K, sqrt(cbind(iso, theo)/pi)  ~ r)
  plot(K, cbind(iso/2, -theo) ~ r)
  plot(K, cbind(iso/2, trans/2) - theo ~ r)

  # test expansion of .x and .y
  plot(K, . ~ .x)
  plot(K, . - theo ~ .x)
  plot(K, .y - theo ~ .x)
  plot(K, sqrt(.y) - sqrt(theo) ~ .x)

  # problems with parsing weird strings in levels(marks(X))
  # noted by Ulf Mehlig

  levels(marks(amacrine)) <- c("Nasticreechia krorluppia", "Homo habilis")
  plot(Kcross(amacrine))
  plot(alltypes(amacrine, "K"))
  plot(alltypes(amacrine, "J"))
  plot(alltypes(amacrine, pcfcross))
})

#
#  Test quirks related to 'alim' attribute

local({
  K <- Kest(cells)
  attr(K, "alim") <- NULL
  plot(K)
  attr(K, "alim") <- c(0, 0.1)
  plot(tail(K))
})

#
# Check that default 'r' vector passes the test for fine spacing

local({
  a <- Fest(cells)
  A <- Fest(cells, r=a$r)
  b <- Hest(heather$coarse)
  B <- Hest(heather$coarse, r=b$r)
  # from Cenk Icos
  X <- runifpoint(100, owin(c(0,3), c(0,10)))
  FX <- Fest(X)
  FXr <- Fest(X, r=FX$r)
  JX <- Jest(X)
})

## various functionality in fv.R

local({
  M <- cbind(1:20, matrix(runif(100), 20, 5))
  A <- as.fv(M)
  fvlabels(A) <- c("r","%s(r)", "%s[A](r)", "%s[B](r)", "%s[C](r)", "%s[D](r)")
  A <- rename.fv(A, "M", quote(M(r)))
  A <- tweak.fv.entry(A, "V1", new.tag="r")
  A[,3] <- NULL
  A$hogwash <- runif(nrow(A))
  #' bind.fv with qualitatively different functions
  GK <- harmonise(G=Gest(cells), K=Kest(cells))
  G <- GK$G
  K <- GK$K
  ss <- c(rep(TRUE, nrow(K)-10), rep(FALSE, 10))
  U <- bind.fv(G, K[ss, ], clip=TRUE)
  #'
  H <- rebadge.as.crossfun(K, "H", "inhom", 1, 2)
  H <- rebadge.as.dotfun(K, "H", "inhom", 3)
})
