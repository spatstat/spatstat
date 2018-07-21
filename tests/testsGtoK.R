##
##    tests/gcc323.R
##
##    $Revision: 1.2 $  $Date: 2015/12/29 08:54:49 $
##
require(spatstat)
local({
  # critical R values that provoke GCC bug #323
  a <- marktable(lansing, R=0.25)
  a <- marktable(lansing, R=0.21)
  a <- marktable(lansing, R=0.20)
  a <- marktable(lansing, R=0.10)
})
#       
#        tests/hobjects.R
#
#   Validity of methods for ppm(... method="ho")
#

require(spatstat)

local({
  set.seed(42)
  fit  <- ppm(cells ~1,         Strauss(0.1), method="ho", nsim=10)
  fitx <- ppm(cells ~offset(x), Strauss(0.1), method="ho", nsim=10)

  a  <- AIC(fit)
  ax <- AIC(fitx)

  f  <- fitted(fit)
  fx <- fitted(fitx)

  p  <- predict(fit)
  px <- predict(fitx)
})


#
# tests/hyperframe.R
#
# test "[.hyperframe" etc
#
#  $Revision: 1.4 $  $Date: 2018/05/15 14:20:38 $
#

require(spatstat)
local({
  lambda <- runif(4, min=50, max=100)
  X <- lapply(as.list(lambda), function(x) { rpoispp(x) })
  h <- hyperframe(lambda=lambda, X=X)
  h$lambda2 <- lambda^2
  h[, "lambda3"] <- lambda^3
  h[, "Y"] <- X
  h[, "X"] <- lapply(X, flipxy)
  h[, c("X", "Y")] <- hyperframe(X=X, Y=X)

  names(h) <- LETTERS[1:5]
  print(h)

  summary(h)
  str(h)
  head(h)
  tail(h)
})


#'     tests/hypotests.R
#'     Hypothesis tests
#' 
#'  $Revision: 1.2 $ $Date: 2018/07/21 03:02:20 $

require(spatstat)
local({
  hopskel.test(redwood, method="MonteCarlo", nsim=5)
  
  berman.test(spiders, "x")
  berman.test(lppm(spiders ~ x), "y")

  #' quadrat test - spatial methods
  a <- quadrat.test(redwood, 3)
  domain(a)
  shift(a, c(1,1))
})
#
#  tests/imageops.R
#
#   $Revision: 1.14 $   $Date: 2018/05/27 05:12:57 $
#

require(spatstat)
local({
  AA <- A <- as.im(owin())
  BB <- B <- as.im(owin(c(1.1, 1.9), c(0,1)))
  Z <- imcov(A, B)
  stopifnot(abs(max(Z) - 0.8) < 0.1)

  Frame(AA) <- Frame(B)
  Frame(BB) <- Frame(A)
  
  ## handling images with 1 row or column
  ycov <- function(x, y) y
  E <- as.im(ycov, owin(), dimyx = c(2,1))
  G <- cut(E, 2)
  H <- as.tess(G)

  E12 <- as.im(ycov, owin(), dimyx = c(1,2))
  G12 <- cut(E12, 2)
  H12 <- as.tess(G12)

  AAA <- as.array(AA)
  EEE <- as.array(E)
  AAD <- as.double(AA)
  EED <- as.double(E)
  aaa <- xtfrm(AAA)
  eee <- xtfrm(E)
  
  ##
  d <- distmap(cells, dimyx=32)
  Z <- connected(d <= 0.06, method="interpreted")

  a <- where.max(d, first=FALSE)
  a <- where.min(d, first=FALSE)

  dx <- raster.x(d)
  dy <- raster.y(d)
  dxy <- raster.xy(d)
  xyZ <- raster.xy(Z, drop=TRUE)

  horosho <- conform.imagelist(cells, list(d, Z))

  #' split.im
  W <- square(1)
  X <- as.im(function(x,y){x}, W)
  Y <- dirichlet(runifpoint(7, W))
  Z <- split(X, as.im(Y))
  
  ## cases of "[.im"
  ee  <- d[simplenet, drop=FALSE]
  eev <- d[simplenet]
  Empty <- cells[FALSE]
  EmptyFun <- ssf(Empty, numeric(0))
  ff <- d[Empty]
  ff <- d[EmptyFun]
  gg <- d[2,]
  gg <- d[,2]
  gg <- d[2:4, 3:5]
  hh <- d[2:4, 3:5, rescue=TRUE]
  if(!is.im(hh)) stop("rectangle was not rescued in [.im")
  ## cases of "[<-.im"
  d[Empty] <- 42
  d[EmptyFun] <- 42
  
  ## smudge() and rasterfilter()
  dd <- smudge(d)

  ## rgb/hsv options
  X <- setcov(owin())
  M <- Window(X)
  Y <- as.im(function(x,y) x, W=M)
  Z <- as.im(function(x,y) y, W=M)
  # convert after rescaling
  RGBscal <- rgbim(X, Y, Z, autoscale=TRUE, maxColorValue=1)
  HSVscal <- hsvim(X, Y, Z, autoscale=TRUE)

  #' miscellaneous
  ZZ <- zapsmall(Z, digits=6)
})
#' indices.R
#' Tests of code for understanding index vectors etc
#' $Revision: 1.1 $ $Date: 2018/03/01 03:38:07 $

require(spatstat)
local({

  a <- grokIndexVector(c(FALSE,TRUE),         10)
  b <- grokIndexVector(rep(c(FALSE,TRUE), 7), 10)
  d <- grokIndexVector(c(2,12),               10)
  e <- grokIndexVector(letters[4:2], nama=letters)
  f <- grokIndexVector(letters[10:1], nama=letters[1:5])
  g <- grokIndexVector(-c(2, 5),              10)
  h <- grokIndexVector(-c(2, 5, 15),          10)

  Nam <- letters[1:10]
  j  <- positiveIndex(-c(2,5), nama=Nam)
  jj <- logicalIndex(-c(2,5), nama=Nam)
  k  <- positiveIndex(-c(2,5), nama=Nam)
  kk <- logicalIndex(-c(2,5), nama=Nam)
  mm <- positiveIndex(c(FALSE,TRUE), nama=Nam)
  nn <- positiveIndex(FALSE, nama=Nam)

  aa <- ppsubset(cells, square(0.1))
})
#'   tests/ippm.R
#'   Tests of 'ippm' class
#'   $Revision: 1.1 $ $Date: 2017/06/06 06:32:00 $

require(spatstat)

local({
  # .......... set up example from help file .................
  nd <- 10
  gamma0 <- 3
  delta0 <- 5
  POW <- 3
  # Terms in intensity
  Z <- function(x,y) { -2*y }
  f <- function(x,y,gamma,delta) { 1 + exp(gamma - delta * x^POW) }
  # True intensity
  lamb <- function(x,y,gamma,delta) { 200 * exp(Z(x,y)) * f(x,y,gamma,delta) }
  # Simulate realisation
  lmax <- max(lamb(0,0,gamma0,delta0), lamb(1,1,gamma0,delta0))
  set.seed(42)
  X <- rpoispp(lamb, lmax=lmax, win=owin(), gamma=gamma0, delta=delta0)
  # Partial derivatives of log f
  DlogfDgamma <- function(x,y, gamma, delta) {
    topbit <- exp(gamma - delta * x^POW)
    topbit/(1 + topbit)
  }
  DlogfDdelta <- function(x,y, gamma, delta) {
    topbit <- exp(gamma - delta * x^POW)
    - (x^POW) * topbit/(1 + topbit)
  }
  # irregular score
  Dlogf <- list(gamma=DlogfDgamma, delta=DlogfDdelta)
  # fit model
  fit <- ippm(X ~Z + offset(log(f)),
              covariates=list(Z=Z, f=f),
              iScore=Dlogf,
              start=list(gamma=1, delta=1),
              nd=nd)

  # ............. test .............................
  Ar <- model.matrix(fit)
  Ai <- model.matrix(fit, irregular=TRUE)
  Zr <- model.images(fit)
  Zi <- model.images(fit, irregular=TRUE)
})#'
#'   tests/Kfuns.R
#'
#'   Various K and L functions and pcf
#'
#'   $Revision: 1.4 $  $Date: 2018/07/13 05:15:00 $
#'

require(spatstat)
myfun <- function(x,y){(x+1) * y }
local({
  #' supporting code
  implemented.for.K(c("border", "bord.modif", "translate", "good", "best"),
                    "mask", TRUE)
  #' Kest special code blocks
  K <- Kest(cells, var.approx=TRUE, ratio=FALSE)
  Z <- distmap(cells) + 1
  Kb <- Kest(cells, correction=c("border","bord.modif"),
             weights=Z, ratio=TRUE)
  Kn <- Kest(cells, correction="none",
             weights=Z, ratio=TRUE)
  #' pcf.ppp special code blocks
  pr  <- pcf(cells, ratio=TRUE, var.approx=TRUE)
  pc  <- pcf(cells, domain=square(0.5))
  pcr <- pcf(cells, domain=square(0.5), ratio=TRUE)
  #' inhomogeneous multitype
  fit <- ppm(amacrine ~ marks)
  K1 <- Kcross.inhom(amacrine, lambdaX=fit)
  K2 <- Kcross.inhom(amacrine, lambdaX=densityfun(amacrine))
  K3 <- Kcross.inhom(amacrine, lambdaX=density(amacrine, at="points"))
  On <- split(amacrine)$on
  Off <- split(amacrine)$off
  K4 <- Kcross.inhom(amacrine, lambdaI=ppm(On), lambdaJ=ppm(Off))
  K5 <- Kcross.inhom(amacrine, correction="bord.modif")
  #' Kmark, markcorr
  X <- runifpoint(100) %mark% runif(100)
  km <- Kmark(X, f=atan2)
  km <- Kmark(X, f1=sin)
  km <- Kmark(X, f="myfun")
  Y <- X %mark% data.frame(u=runif(100), v=runif(100))
  mk <- markcorr(Y)
  #'
  rr <- rep(0.1, npoints(cells))
  eC <- edge.Ripley(cells, rr)
  eI <- edge.Ripley(cells, rr, method="interpreted")
  if(max(abs(eC-eI)) > 0.1)
    stop("Ripley edge correction results do not match")

  a <- rmax.Ripley(square(1))
  a <- rmax.Ripley(as.polygonal(square(1)))
  a <- rmax.Ripley(letterR)
})
#
# tests/kppm.R
#
# $Revision: 1.23 $ $Date: 2018/07/21 00:49:50 $
#
# Test functionality of kppm that depends on RandomFields
# Test update.kppm for old style kppm objects

require(spatstat)
local({

 fit <- kppm(redwood ~1, "Thomas") # sic
 fitx <- update(fit, ~ . + x)
 fitM <- update(fit, clusters="MatClust")
 fitC <- update(fit, cells)
 fitCx <- update(fit, cells ~ x)

 #'
 Wsub <- owin(c(0, 0.5), c(-0.5, 0))
 fitsub <- kppm(redwood ~1, "Thomas", subset=Wsub)
 fitsub
 
 #' various methods
 ff <- as.fv(fitx)
 Y <- simulate(fitx, seed=42)[[1]]
 uu <- unitname(fitx)
 unitname(fitCx) <- "furlong"
 mo <- model.images(fitCx)
 
 # vcov.kppm different algorithms
 vc  <- vcov(fitx)
 vc2 <- vcov(fitx, fast=TRUE)
 vc3 <- vcov(fitx, fast=TRUE, splitup=TRUE)
 vc4 <- vcov(fitx,            splitup=TRUE)

 ## other code blocks
 a <- varcount(fitx, function(x,y){x+1}) # always positive
 a <- varcount(fitx, function(x,y){y-1}) # always negative
 a <- varcount(fitx, function(x,y){x+y}) # positive or negative
 
 # improve.kppm
 fitI <- update(fit, improve.type="quasi")
 fitxI <- update(fitx, improve.type="quasi")
 # vcov.kppm
 vcI <- vcov(fitxI)
 
  # plot.kppm including predict.kppm
 fitMC <- kppm(redwood ~ x, "Thomas")
 fitCL <- kppm(redwood ~ x, "Thomas", method="c")
 fitPA <- kppm(redwood ~ x, "Thomas", method="p")
 plot(fitMC)
 plot(fitCL)
 plot(fitPA)

 # fit with composite likelihood method [thanks to Abdollah Jalilian]
 fut <- kppm(redwood ~ x, "VarGamma", method="clik2", nu.ker=-3/8)
 kfut <- as.fv(fut)
 
 if(require(RandomFields)) {
   fit0 <- kppm(redwood ~1, "LGCP")
   is.poisson(fit0)
   Y0 <- simulate(fit0)[[1]]
   stopifnot(is.ppp(Y0))

   ## fit LGCP using K function: slow
   fit1 <- kppm(redwood ~x, "LGCP",
                covmodel=list(model="matern", nu=0.3),
                control=list(maxit=3))
   Y1 <- simulate(fit1)[[1]]
   stopifnot(is.ppp(Y1))

   ## fit LGCP using pcf
   fit1p <- kppm(redwood ~x, "LGCP",
                 covmodel=list(model="matern", nu=0.3),
                 statistic="pcf")
   Y1p <- simulate(fit1p)[[1]]
   stopifnot(is.ppp(Y1p))

   ## .. and using different fitting methods
   fit1pClik <- update(fit1p, method="clik")
   fit1pPalm <- update(fit1p, method="palm")
   
   ## image covariate (a different code block) 
   xx <- as.im(function(x,y) x, Window(redwood))
   fit1xx <- update(fit1p, . ~ xx, data=solist(xx=xx))
   Y1xx <- simulate(fit1xx)[[1]]
   stopifnot(is.ppp(Y1xx))
   fit1xxVG <- update(fit1xx, clusters="VarGamma", nu=-1/4)
   Y1xxVG <- simulate(fit1xxVG)[[1]]
   stopifnot(is.ppp(Y1xxVG))
   
   # ... and Abdollah's code

   fit2 <- kppm(redwood ~x, cluster="Cauchy", statistic="K")
   Y2 <- simulate(fit2)[[1]]
   stopifnot(is.ppp(Y2))

 }

})

local({
  #' various code blocks
  fut <- kppm(redwood, ~x)
  fet <- update(fut, redwood3)
  fot <- update(fut, trend=~y)
  fit <- kppm(redwoodfull ~ x)
  Y <- simulate(fit, window=redwoodfull.extra$regionII)
  gut <- improve.kppm(fit, type="wclik1")
  gut <- improve.kppm(fit, vcov=TRUE, fast.vcov=TRUE, save.internals=TRUE)
  hut <- kppm(redwood ~ x, method="clik", weightfun=NULL)
  hut <- kppm(redwood ~ x, method="palm", weightfun=NULL)
})

local({
  K <- Kest(redwood)
  a <- matclust.estK(K)
  a <- thomas.estK(K)
  a <- cauchy.estK(K)
  a <- lgcp.estK(K)
  g <- pcf(redwood)
  a <- matclust.estpcf(g)
  a <- thomas.estpcf(g)
  a <- cauchy.estpcf(g)
  a <- lgcp.estpcf(g)
})

local({
  #'  experimental
  spatstat.options(kppm.canonical=TRUE, kppm.adjusted=TRUE)
  futTT1 <- kppm(redwood)
  futTT2 <- kppm(redwood, method="palm")
  futTT3 <- kppm(redwood, method="clik2")
  spatstat.options(kppm.canonical=TRUE, kppm.adjusted=FALSE)
  futTF1 <- kppm(redwood)
  futTF2 <- kppm(redwood, method="palm")
  futTF3 <- kppm(redwood, method="clik2")
  spatstat.options(kppm.canonical=FALSE, kppm.adjusted=TRUE)
  futFT1 <- kppm(redwood)
  futFT2 <- kppm(redwood, method="palm")
  futFT3 <- kppm(redwood, method="clik2")
  spatstat.options(kppm.canonical=FALSE, kppm.adjusted=FALSE)
  futFF1 <- kppm(redwood)
  futFF2 <- kppm(redwood, method="palm")
  futFF3 <- kppm(redwood, method="clik2")
})

reset.spatstat.options()
  
