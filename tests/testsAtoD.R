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
#'
#'   tests/cluck.R
#'
#'   Tests of "click*" functions
#'   using queueing feature of spatstatLocator
#'
#'   $Revision: 1.1 $ $Date: 2018/10/19 10:03:23 $

require(spatstat)
local({
  #' clickppp
  spatstat.utils::queueSpatstatLocator(runif(5), runif(5))
  XA <- clickppp(hook=square(0.5))
  spatstat.utils::queueSpatstatLocator(runif(6), runif(6))
  XB <- clickppp(n=3, types=c("a", "b"))
  #' clickbox
  spatstat.utils::queueSpatstatLocator(runif(2), runif(2))
  BB <- clickbox()
  #' clickdist
  spatstat.utils::queueSpatstatLocator(runif(2), runif(2))
  dd <- clickdist()
  #' clickpoly
  hex <- vertices(disc(radius=0.4, centre=c(0.5, 0.5), npoly=6))
  spatstat.utils::queueSpatstatLocator(hex)
  PA <- clickpoly()
  holy <- vertices(disc(radius=0.2, centre=c(0.5, 0.5), npoly=6))
  holy <- lapply(holy, rev)
  spatstat.utils::queueSpatstatLocator(concatxy(hex, holy))
  PB <- clickpoly(np=2, nv=6)
  #' clicklpp
  Y <- coords(runiflpp(6, simplenet))
  spatstat.utils::queueSpatstatLocator(Y)
  XL <- clicklpp(simplenet)
  spatstat.utils::queueSpatstatLocator(Y)
  XM <- clicklpp(simplenet, n=3, types=c("a", "b"))
})
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
# $Revision: 1.6 $ $Date: 2018/10/07 09:58:42 $

require(spatstat)

local({
  eps <- .Machine$double.eps * 4

  checkagree <- function(A, B, blurb) {
    maxerr <- max(abs(A-B))
    cat("Discrepancy", maxerr, "for", blurb, fill=TRUE)
    if(maxerr > eps) 
      stop(paste("Algorithms for", blurb, "disagree"))
    return(TRUE)
  }

  ## pairdist.ppp
  set.seed(190901)
  X <- rpoispp(42)
  dC <- pairdist(X, method="C")
  dR <- pairdist(X, method="interpreted")
  checkagree(dC, dR, "pairdist()")

  dCp <- pairdist(X, periodic=TRUE, method="C")
  dRp <- pairdist(X, periodic=TRUE, method="interpreted")
  checkagree(dCp, dRp, "pairdist(periodic=TRUE)")

  dCp2 <- pairdist(X, periodic=TRUE, squared=TRUE, method="C")
  dRp2 <- pairdist(X, periodic=TRUE, squared=TRUE, method="interpreted")
  checkagree(dCp2, dRp2, "pairdist(periodic=TRUE, squared=TRUE)")

  ## crossdist.ppp
  Y <- rpoispp(42)
  dC <- crossdist(X, Y, method="C")
  dR <- crossdist(X, Y, method="interpreted")
  checkagree(dC, dR, "crossdist()")

  dC <- crossdist(X, Y, periodic=TRUE, method="C")
  dR <- crossdist(X, Y, periodic=TRUE, method="interpreted")
  checkagree(dC, dR, "crossdist(periodic=TRUE)")

  dC2 <- crossdist(X, Y, periodic=TRUE, squared=TRUE, method="C")
  dR2 <- crossdist(X, Y, periodic=TRUE, squared=TRUE, method="interpreted")
  checkagree(dC2, dR2, "crossdist(periodic=TRUE, squared=TRUE)")

  # nndist.ppp
  nnC <- nndist(X, method="C")
  nnI <- nndist(X, method="interpreted")
  checkagree(nnC, nnI, "nndist()")

  nn3C <- nndist(X, k=3, method="C")
  nn3I <- nndist(X, k=3, method="interpreted")
  checkagree(nn3C, nn3I, "nndist(k=3)")

  # nnwhich.ppp
  nwC <- nnwhich(X, method="C")
  nwI <- nnwhich(X, method="interpreted")
  checkagree(nwC, nwI, "nnwhich()")

  nw3C <- nnwhich(X, k=3, method="C")
  nw3I <- nnwhich(X, k=3, method="interpreted")
  checkagree(nw3C, nw3I, "nnwhich(k=3)")

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
#'  $Revision: 1.33 $  $Date: 2018/10/21 10:09:11 $
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

  ## Kinhom 
  lam <- density(redwood)
  K <- Kinhom(redwood, lam)
  
  lamX <- density(redwood, at="points")
  KX <- Kinhom(redwood, lamX)

  ## test all code cases of new 'relrisk.ppp' algorithm
  pants <- function(..., X=ants, sigma=100) {
    a <- relrisk(X, sigma=sigma, se=TRUE, ...)
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

  pants(sigma=Inf)
  
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
  stroke <- function(..., Y = longleaf) {
    Z <- Smooth(Y, ..., at="pixels")
    Z <- Smooth(Y, ..., at="points")
    f <- Smoothfun(Y, ...)
    f(120, 80)
    f(Y[1:2])
    f(Y[FALSE])
    return(invisible(NULL))
  }
  stroke()
  stroke(5, diggle=TRUE)
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
    Z <- Smooth(Y, ..., at="points")
    f <- Smoothfun(Y, ...)
    f(4, 1)
    f(Y[1:2])
    f(Y[FALSE])
    return(invisible(NULL))
  }
  strike()
  strike(varcov=diag(c(1.2, 2.1)))
  strike(1.5, weights=runif(npoints(finpines)))
  strike(1.5, weights=expression(y))
  strike(1e-6)

  strike(sigma=Inf)

  ## detect special cases
  Smooth(longleaf[FALSE])
  Xconst <- cells %mark% 1
  Smooth(Xconst, 0.1)
  Smooth(Xconst, 0.1, at="points")
  Smooth(cells %mark% runif(42), sigma=Inf)
  Smooth(cells %mark% runif(42), sigma=Inf, at="points")
  Smooth(cells %mark% runif(42), sigma=Inf, at="points", leaveoneout=FALSE)
  
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

local({
  ## idw
  Z <- idw(longleaf, power=4)
  ZX <- idw(longleaf, power=4, at="points")
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
#'   $Revision: 1.2 $  $Date: 2018/10/26 08:12:26 $

require(spatstat)
local({
  X <- runifpoint(10)
  Y <- runifpoint(10)
  d <- pppdist(X, Y, type="ace", show.rprimal=TRUE)
  d2 <- pppdist(X, Y, type="spa", ccode=FALSE)
  d3 <- pppdist(X, Y, type="mat", ccode=TRUE, auction=FALSE)
  m <- pppdist.mat(X, Y, q=Inf, cutoff=0.001)
  a <- matchingdist(d, type="ace")
  b <- matchingdist(d, type="mat")
})



#'
#'   tests/duplicity.R
#'
#'  Tests of duplicated/multiplicity code
#'
#' $Revision: 1.2 $ $Date: 2018/11/02 01:23:09 $

require(spatstat)
local({
   X <- ppp(c(1,1,0.5,1), c(2,2,1,2), window=square(3), check=FALSE)
   m <- multiplicity(X)
   Y <- X %mark% factor(letters[c(3,2,4,3)])
   mf <- multiplicity(Y)
   Z <- Y %mark% matrix(c(3,2,4,3), 4, 2)
   mm <- multiplicity(Z)
   marks(Z) <- as.data.frame(marks(Z))
   mz <- multiplicity(Z)
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
})
