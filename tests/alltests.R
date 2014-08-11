# nndist.R
# Check that nndist and nnwhich give
# results consistent with direct calculation from pairdist

# Similarly for nncross and distfun

require(spatstat)
local({
  eps <- sqrt(.Machine$double.eps)
  f <- function(mat,k) { apply(mat, 1, function(z,n) { sort(z)[n]  }, n=k+1) }
  g <- function(mat,k) { apply(mat, 1, function(z,n) { order(z)[n] }, n=k+1) }

  # Two dimensions

  X <- runifpoint(42)

  nn <- nndist(X)
  nnP <- f(pairdist(X), 1)
  if(any(abs(nn - nnP) > eps))
    stop("nndist.ppp does not agree with pairdist")

  nn5 <- nndist(X, k=5)
  nn5P <- f(pairdist(X), 5)
  if(any(abs(nn5 - nn5P) > eps))
    stop("nndist.ppp(k=5) does not agree with pairdist")

  nw <- nnwhich(X)
  nwP <- g(pairdist(X), 1)
  if(any(nw != nwP))
    stop("nnwhich.ppp does not agree with pairdist")

  nw5 <- nnwhich(X, k=5)
  nw5P <- g(pairdist(X), 5)
  if(any(nw5 != nw5P))
    stop("nnwhich.ppp(k=5) does not agree with pairdist")

  # Three dimensions

  X <- runifpoint3(42)

  nn <- nndist(X)
  nnP <- f(pairdist(X), 1)
  if(any(abs(nn - nnP) > eps))
    stop("nndist.pp3 does not agree with pairdist")

  nn5 <- nndist(X, k=5)
  nn5P <- f(pairdist(X), 5)
  if(any(abs(nn5 - nn5P) > eps))
    stop("nndist.pp3(k=5) does not agree with pairdist")

  nw <- nnwhich(X)
  nwP <- g(pairdist(X), 1)
  if(any(nw != nwP))
    stop("nnwhich.pp3 does not agree with pairdist")

  nw5 <- nnwhich(X, k=5)
  nw5P <- g(pairdist(X), 5)
  if(any(nw5 != nw5P))
    stop("nnwhich.pp3(k=5) does not agree with pairdist")

  # m dimensions

  X <- runifpointx(42, boxx(c(0,1),c(0,1),c(0,1),c(0,1)))

  nn <- nndist(X)
  nnP <- f(pairdist(X), 1)
  if(any(abs(nn - nnP) > eps))
    stop("nndist.ppx does not agree with pairdist")

  nn5 <- nndist(X, k=5)
  nn5P <- f(pairdist(X), 5)
  if(any(abs(nn5 - nn5P) > eps))
    stop("nndist.ppx(k=5) does not agree with pairdist")
  
  nw <- nnwhich(X)
  nwP <- g(pairdist(X), 1)
  if(any(nw != nwP))
    stop("nnwhich.ppx does not agree with pairdist")

  nw5 <- nnwhich(X, k=5)
  nw5P <- g(pairdist(X), 5)
  if(any(nw5 != nw5P))
    stop("nnwhich.ppx(k=5) does not agree with pairdist")

  #### nncross in two dimensions
  X <- runifpoint(42)
  Y <- runifpoint(42, win=owin(c(1,2),c(1,2)))

  # default nncross
  nc <- nncross(X,Y)
  ncd <- nc$dist
  ncw <- nc$which
  cd <- crossdist(X,Y)
  cdd <- apply(cd, 1, min)
  cdw <- apply(cd, 1, which.min)
  if(any(abs(ncd - cdd) > eps))
    stop("nncross()$dist does not agree with apply(crossdist(), 1, min)")
  if(any(ncw != cdw))
    stop("nncross()$which does not agree with apply(crossdist(), 1, which.min)")

  # sort on x
  nc <- nncross(X,Y, sortby="x")
  ncd <- nc$dist
  ncw <- nc$which
  if(any(abs(ncd - cdd) > eps))
    stop("nncross(sortby=x)$dist does not agree with apply(crossdist(), 1, min)")
  if(any(ncw != cdw))
    stop("nncross(sortby=x)$which does not agree with apply(crossdist(), 1, which.min)")

  # pre-sorted on x
  Y <- Y[order(Y$x)]
  nc <- nncross(X,Y, is.sorted.Y=TRUE, sortby="x")
  ncd <- nc$dist
  ncw <- nc$which
  cd <- crossdist(X,Y)
  cdd <- apply(cd, 1, min)
  cdw <- apply(cd, 1, which.min)
  if(any(abs(ncd - cdd) > eps))
    stop("For sorted data, nncross()$dist does not agree with apply(crossdist(), 1, min)")
  if(any(ncw != cdw))
    stop("For sorted data, nncross()$which does not agree with apply(crossdist(), 1, which.min)")

  # sanity check for nncross with k > 1
  ndw <- nncross(X, Y, k=1:4)$which
  if(any(is.na(ndw)))
    stop("NA's returned by nncross.ppp(k > 1)$which")
  ndw <- nncross(X, Y, k=1:4, what="which")
  if(any(is.na(ndw)))
    stop("NA's returned by nncross.ppp(k > 1, what='which')")
  
  # test of correctness for nncross with k > 1
  flipcells <- flipxy(cells)
  calcwhich <- nncross(cells, flipcells, k=1:4, what="which")
  truewhich <- t(apply(crossdist(cells,flipcells), 1, order))[,1:4]
  if(any(calcwhich != truewhich))
    stop("nncross(k > 1) gives wrong answer")
  
  # test of agreement between nngrid.h and knngrid.h
  #    dimyx=23 (found by trial-and-error) ensures that there are no ties 
  a <- as.matrix(nnmap(cells, what="which", dimyx=23))
  b <- as.matrix(nnmap(cells, what="which", dimyx=23, k=1:2)[[1]])
  if(any(a != b))
    stop("algorithms in nngrid.h and knngrid.h disagree")
})


require(spatstat)

local({
  Y <- split(urkiola)
  B <- Y$birch
  O <- Y$oak
  B.lam <- predict (ppm(B, ~polynom(x,y,2)), type="trend")
  O.lam <- predict (ppm(O, ~polynom(x,y,2)), type="trend")

  Kinhom(B, lambda=B.lam, correction="iso")
  Kinhom(B, lambda=B.lam, correction="border")

  Kcross.inhom(urkiola, i="birch", j="oak", B.lam, O.lam)
  Kcross.inhom(urkiola, i="birch", j="oak", B.lam, O.lam, correction = "iso")
  Kcross.inhom(urkiola, i="birch", j="oak", B.lam, O.lam, correction = "border")
})


require(spatstat)
local({
  # critical R values that provoke GCC bug #323
  a <- marktable(lansing, R=0.25)
  a <- marktable(lansing, R=0.21)
  a <- marktable(lansing, R=0.20)
  a <- marktable(lansing, R=0.10)
})
# tests/NAinCov.R
# Testing the response to the presence of NA's in covariates

require(spatstat)
local({
  X <- runifpoint(42)
  Y <- as.im(function(x,y) { x+y }, owin())
  Y[owin(c(0.2,0.4),c(0.2,0.4))] <- NA
  # fit model: should produce a warning but no failure
  misfit <- ppm(X, ~Y, covariates=list(Y=Y))
  # prediction 
  Z <- predict(misfit, type="trend")
  Z <- predict(misfit, type="se")
  # covariance matrix: all should be silent
  v <- vcov(misfit)
  ss <- vcov(misfit, what="internals")
  NULL
})





# tests for agreement between C and interpreted code
# for interpoint distances

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
})
# ppmBadData.R
# $Revision: 1.4 $ $Date: 2011/12/05 07:29:16 $

# Testing robustness of ppm and support functions
# when data are rubbish

require(spatstat)
local({
# ---------------------------------------------------
# from Rolf: very large proportion of data is NA
  SEED <- 42
  K <- 101
  A <- 500
  X <- seq(0, A, length=K)
  G <- expand.grid(x=X, y=X)
  FOO <- function(x,y) { sin(x)^2 + cos(y)^2 }
  M1 <- im(matrix(FOO(G$x, G$y), K, K), xcol=X, yrow=X)
  M <- im(matrix(FOO(G$x, G$y), K, K))
  BAR <- function(x) { exp(-6.618913 + 5.855337 * x - 8.432483 * x^2) }
  V <- im(BAR(M$v), xcol=X, yrow=X)
  # V <- eval.im(exp(-6.618913 + 5.855337 * M - 8.432483 * M^2))
  set.seed(SEED)
  Y <- rpoispp(V)
  fY <- ppm(Y, ~cv + I(cv^2), covariates=list(cv=M), correction="translate")
  diagnose.ppm(fY)
  lurking(fY, covariate=as.im(function(x,y){x}, square(A)), type="raw")
})

# --------------------------------------------------------
# from Andrew Bevan: numerical overflow, ill-conditioned Fisher information

local({
  SEED <- 42

  nongranite<- owin(poly = list(x = c(0, 8500, 7000, 6400, 6400, 6700, 7000, 7200, 7300, 8000, 8100, 8800, 9500, 10000, 10000, 0), y = c(0, 0, 2000, 3800, 4000, 5000, 6500, 7400, 7500, 8000, 8100, 9000, 9500, 9600, 10000, 10000)))

  #Trend on raster grid
  rain <- as.im(X=function(x,y) { x^2 + y^2 }, W=nongranite, dimyx=100)

  #Generate a point pattern via a Lennard-Jones process
  set.seed(SEED)
  mod4<- rmhmodel(cif="lennard",
                par=list(beta=1, sigma=250, epsilon=2.2),
                trend=rain, w=nongranite)
  ljtr<- rmh(mod4, start=list(n.start=80), control=list(p=1, nrep=1e5))

  #Fit a point process model to the pattern with rain as a covariate
  # NOTE INCORRECT TREND FORMULA
  ljtrmod <- ppm(ljtr, trend= ~ Z, interaction=NULL, covariates=list(Z=rain))
  ss <- summary(ljtrmod)
})

local({
  # From Ege
  # Degenerate but non-null argument 'covariates'
  xx <- list()
  names(xx) <- character(0)
  fit <- ppm(cells, ~x, covariates = xx)
  st <- summary(fit) 
})
#
#   tests/ppmgam.R
#
#   Test ppm with use.gam=TRUE
#
#   $Revision: 1.2 $  $Date: 2013/03/01 08:46:34 $
#

require(spatstat)
local({
  fit <- ppm(nztrees, ~s(x,y), use.gam=TRUE)
  mm <- model.matrix(fit)
  mf <- model.frame(fit)
  v <- vcov(fit)
})

#
#  ppmlogi.R
#
# Tests of ppm(method='logi')
#
# $Revision: 1.2 $  Date$
#

require(spatstat)
local({
  fit <- ppm(cells, ~x, method="logi")
  f <- fitted(fit)
  p <- predict(fit)
  fitS <- ppm(cells, ~x, Strauss(0.08), method="logi")
  fS <- fitted(fitS)
  pS <- predict(fitS)
  if(spatstat.options("allow.logi.influence")) {
    a <- leverage(fit)
    b <- influence(fit)
    d <- dfbetas(fit)
    aS <- leverage(fitS)
    bS <- influence(fitS)
    dS <- dfbetas(fitS)
  }
})
# ppmmarkorder.R
# $Revision: 1.2 $  $Date: 2011/12/05 07:29:16 $
# Test that predict.ppm, plot.ppm and plot.fitin
# tolerate marks with levels that are not in alpha order
#
require(spatstat)
local({
  X <- amacrine
  levels(marks(X)) <- c("ZZZ", "AAA")
  fit <- ppm(X, ~marks, MultiStrauss(c("ZZZ","AAA"), matrix(0.06, 2, 2)))
  aa <- predict(fit, type="trend")
  bb <- predict(fit, type="cif")
  plot(fit)
  plot(fitin(fit))
})


#
#   tests/ppmscope.R
#
#   Test things that might corrupt the internal format of ppm objects
#
#   $Revision: 1.3 $  $Date: 2011/12/05 07:29:16 $
#
#   (1) Scoping problem that can arise when ppm splits the data

require(spatstat)
local({
  fit <- ppm(bei, ~elev, covariates=bei.extra)
  mm <- model.matrix(fit)

  #   (2) Fast update mechanism

  fit1 <- ppm(cells, ~x+y, Strauss(0.07))
  fit2 <- update(fit1, ~y)
  fit3 <- update(fit2, ~x)
})

#
#   tests/ppmtricks.R
#
#   Test backdoor exits and hidden options in ppm
#
#   $Revision: 1.1 $  $Date: 2013/06/17 06:54:51 $
#
require(spatstat)
local({

  # (1) skip.border
  
  fit <- ppm(cells, ~1, Strauss(0.1), skip.border=TRUE)

})

#
# tests/prediction.R
#
# Things that might go wrong with predict()
#
#  $Revision: 1.1 $ $Date: 2013/11/12 16:06:11 $
#

require(spatstat)

local({
  # test of 'covfunargs'
  f <- function(x,y,a){ y - a }
  fit <- ppm(cells, ~x + f, covariates=list(f=f), covfunargs=list(a=1/2))
  p <- predict(fit)
})

#
#  tests/rmhAux.R
#
#  $Revision: 1.1 $  $Date: 2013/02/18 10:41:27 $
#
#  For interactions which maintain 'auxiliary data',
#  verify that the auxiliary data are correctly updated.
#
#  To do this we run rmh with nsave=1 so that the point pattern state
#  is saved after every iteration, then the algorithm is restarted,
#  and the auxiliary data are re-initialised. The final state must agree with
#  the result of simulation without saving.
# ----------------------------------------------------

require(spatstat)

local({

   # Geyer:
   mod <- list(cif="geyer",
               par=list(beta=1.25,gamma=1.6,r=0.2,sat=4.5),
               w=square(10))

   set.seed(42)
   X.nosave <- rmh(model=mod,
                   start=list(n.start=50),
                   control=list(nrep=1e3, periodic=FALSE, expand=1))
   set.seed(42)
   X.save <- rmh(model=mod,
                 start=list(n.start=50),
                 control=list(nrep=1e3, periodic=FALSE, expand=1,
                   nburn=0, nsave=1))

   stopifnot(npoints(X.save) == npoints(X.nosave))
   stopifnot(max(nncross(X.save, X.nosave)$dist) == 0)
   stopifnot(max(nncross(X.nosave, X.save)$dist) == 0)
})
# Test examples for rmh.default
# run to reasonable length
# and with tests for validity added
# ----------------------------------------------------

require(spatstat)

local({
if(!exists("nr"))
   nr   <- 5e3

spatstat.options(expand=1.1)
   
   # Strauss process.
   mod01 <- list(cif="strauss",par=list(beta=2,gamma=0.2,r=0.7),
                 w=c(0,10,0,10))
   X1.strauss <- rmh(model=mod01,start=list(n.start=80),
                     control=list(nrep=nr))

   # Strauss process, conditioning on n = 80:
   X2.strauss <- rmh(model=mod01,start=list(n.start=80),
                     control=list(p=1,nrep=nr))
   stopifnot(X2.strauss$n == 80)

   # test tracking mechanism
   X1.strauss <- rmh(model=mod01,start=list(n.start=80),
                     control=list(nrep=nr), track=TRUE)
   X2.strauss <- rmh(model=mod01,start=list(n.start=80),
                         control=list(p=1,nrep=nr), track=TRUE)
   
   # Hard core process:
   mod02 <- list(cif="hardcore",par=list(beta=2,hc=0.7),w=c(0,10,0,10))
   X3.hardcore <- rmh(model=mod02,start=list(n.start=60),
                     control=list(nrep=nr))
   
   # Strauss process equal to pure hardcore:
   mod02 <- list(cif="strauss",par=list(beta=2,gamma=0,r=0.7),w=c(0,10,0,10))
   X3.strauss <- rmh(model=mod02,start=list(n.start=60),
                     control=list(nrep=nr))
   
   # Strauss process in a polygonal window.
   x     <- c(0.55,0.68,0.75,0.58,0.39,0.37,0.19,0.26,0.42)
   y     <- c(0.20,0.27,0.68,0.99,0.80,0.61,0.45,0.28,0.33)
   mod03 <- list(cif="strauss",par=list(beta=2000,gamma=0.6,r=0.07),
                w=owin(poly=list(x=x,y=y)))
   X4.strauss <- rmh(model=mod03,start=list(n.start=90),
                     control=list(nrep=nr))
   
   # Strauss process in a polygonal window, conditioning on n = 42.
   X5.strauss <- rmh(model=mod03,start=list(n.start=42),
                     control=list(p=1,nrep=nr))
   stopifnot(X5.strauss$n == 42)

   # Strauss process, starting off from X4.strauss, but with the
   # polygonal window replace by a rectangular one.  At the end,
   # the generated pattern is clipped to the original polygonal window.
   xxx <- X4.strauss
   xxx$window <- as.owin(c(0,1,0,1))
   X6.strauss <- rmh(model=mod03,start=list(x.start=xxx),
                     control=list(nrep=nr))
   
   # Strauss with hardcore:
   mod04 <- list(cif="straush",par=list(beta=2,gamma=0.2,r=0.7,hc=0.3),
                w=c(0,10,0,10))
   X1.straush <- rmh(model=mod04,start=list(n.start=70),
                     control=list(nrep=nr))
   
   # Another Strauss with hardcore (with a perhaps surprising result):
   mod05 <- list(cif="straush",par=list(beta=80,gamma=0.36,r=45,hc=2.5),
                w=c(0,250,0,250))
   X2.straush <- rmh(model=mod05,start=list(n.start=250),
                     control=list(nrep=nr))
   
   # Pure hardcore (identical to X3.strauss).
   mod06 <- list(cif="straush",par=list(beta=2,gamma=1,r=1,hc=0.7),
                w=c(0,10,0,10))
   X3.straush <- rmh(model=mod06,start=list(n.start=60),
                     control=list(nrep=nr))

   # Area-interaction, inhibitory
   mod.area <- list(cif="areaint",par=list(beta=2,eta=0.5,r=0.5), w=square(10))
   X.area <- rmh(model=mod.area,start=list(n.start=60),
                 control=list(nrep=nr))

   # Area-interaction, clustered
   mod.area2 <- list(cif="areaint",par=list(beta=2,eta=1.5,r=0.5), w=square(10))
   X.area2 <- rmh(model=mod.area2,start=list(n.start=60),
                 control=list(nrep=nr))

   # Area-interaction close to hard core
   set.seed(42)
   mod.area0 <- list(cif="areaint",par=list(beta=2,eta=1e-300,r=0.35),
                     w=square(10))
   X.area0 <- rmh(model=mod.area0,start=list(x.start=X3.hardcore),
                 control=list(nrep=nr))
   stopifnot(nndist(X.area0) > 0.6)
   
   # Soft core:
   w    <- c(0,10,0,10)
   mod07 <- list(cif="sftcr",par=list(beta=0.8,sigma=0.1,kappa=0.5),
                w=c(0,10,0,10))
   X.sftcr <- rmh(model=mod07,start=list(n.start=70),
                  control=list(nrep=nr))
   
   # Diggle, Gates, and Stibbard:
   mod12 <- list(cif="dgs",par=list(beta=3600,rho=0.08),w=c(0,1,0,1))
   X.dgs <- rmh(model=mod12,start=list(n.start=300),
                control=list(nrep=nr))
   
   # Diggle-Gratton:
   mod13 <- list(cif="diggra",
                 par=list(beta=1800,kappa=3,delta=0.02,rho=0.04),
                 w=square(1))
   X.diggra <- rmh(model=mod13,start=list(n.start=300),
                   control=list(nrep=nr))
   
   # Geyer:
   mod14 <- list(cif="geyer",par=list(beta=1.25,gamma=1.6,r=0.2,sat=4.5),
                 w=c(0,10,0,10))
   X1.geyer <- rmh(model=mod14,start=list(n.start=200),
                   control=list(nrep=nr))
   
   # Geyer; same as a Strauss process with parameters
   # (beta=2.25,gamma=0.16,r=0.7):
   
   mod15 <- list(cif="geyer",par=list(beta=2.25,gamma=0.4,r=0.7,sat=10000),
                 w=c(0,10,0,10))
   X2.geyer <- rmh(model=mod15,start=list(n.start=200),
                   control=list(nrep=nr))
   
   mod16 <- list(cif="geyer",par=list(beta=8.1,gamma=2.2,r=0.08,sat=3))
   data(redwood)
   X3.geyer <- rmh(model=mod16,start=list(x.start=redwood),
                   control=list(periodic=TRUE,nrep=nr))
   
   # Geyer, starting from the redwood data set, simulating
   # on a torus, and conditioning on n:
   X4.geyer <- rmh(model=mod16,start=list(x.start=redwood),
                   control=list(p=1,periodic=TRUE,nrep=nr))

   # Lookup (interaction function h_2 from page 76, Diggle (2003)):
      r <- seq(from=0,to=0.2,length=101)[-1] # Drop 0.
      h <- 20*(r-0.05)
      h[r<0.05] <- 0
      h[r>0.10] <- 1
      mod17 <- list(cif="lookup",par=list(beta=4000,h=h,r=r),w=c(0,1,0,1))
      X.lookup <- rmh(model=mod17,start=list(n.start=100),
                      control=list(nrep=nr))
                   
   # Strauss with trend
   tr <- function(x,y){x <- x/250; y <- y/250;
   			   exp((6*x + 5*y - 18*x^2 + 12*x*y - 9*y^2)/6)
                         }
   beta <- 0.3
   gmma <- 0.5
   r    <- 45
   tr3   <- function(x,y){x <- x/250; y <- y/250;
   			   exp((6*x + 5*y - 18*x^2 + 12*x*y - 9*y^2)/6)
                         }
                         # log quadratic trend
   mod17 <- list(cif="strauss",par=list(beta=beta,gamma=gmma,r=r),w=c(0,250,0,250),
                 trend=tr3)
   X1.strauss.trend <- rmh(model=mod17,start=list(n.start=90),
                           control=list(nrep=nr))

})
# Things which should cause an error
require(spatstat)

local({
if(!exists("nv"))
  nv <- 0
if(!exists("nr"))
  nr   <- 1e3

# Strauss with zero intensity and p = 1
mod0S <- list(cif="strauss",par=list(beta=0,gamma=0.6,r=0.7), w = square(3))
out <- try(X0S   <- rmh(model=mod0S,start=list(n.start=80),
               control=list(p=1,nrep=nr,nverb=nv),verbose=FALSE))
if(!inherits(out, "try-error"))
  stop("Error not trapped (Strauss with zero intensity and p = 1) in tests/rmhErrors.R")
})


#
# tests/rmhExpand.R
#
# test decisions about expansion of simulation window
#
#  $Revision: 1.2 $  $Date: 2011/12/05 07:29:16 $
#

require(spatstat)
local({
fit <- ppm(cells, ~x)

# check rmhmodel.ppm
mod <- rmhmodel(fit)
wsim <- as.rectangle(mod$trend)
if(!identical(wsim, as.owin(cells)))
  stop("Expansion occurred improperly in rmhmodel.ppm")
})


#
#  tests of rmh, running multitype point processes
#
require(spatstat)

local({
if(!exists("nr"))
   nr   <- 5e3

if(!exists("nv"))
   nv   <- 0

spatstat.options(expand=1.1)

   # Multitype Poisson
   modp2 <- list(cif="poisson",
                 par=list(beta=2), types=letters[1:3], w = square(10))
   Xp2 <- rmh(modp2, start=list(n.start=0), control=list(p=1))
    
   # Multitype Strauss:
   beta <- c(0.027,0.008)
   gmma <- matrix(c(0.43,0.98,0.98,0.36),2,2)
   r    <- matrix(c(45,45,45,45),2,2)
   mod08 <- list(cif="straussm",par=list(beta=beta,gamma=gmma,radii=r),
                w=c(0,250,0,250))
   X1.straussm <- rmh(model=mod08,start=list(n.start=80),
                      control=list(ptypes=c(0.75,0.25),nrep=nr,nverb=nv))
   

   # Multitype Strauss conditioning upon the total number
   # of points being 80:
   X2.straussm <- rmh(model=mod08,start=list(n.start=80),
                      control=list(p=1,ptypes=c(0.75,0.25),nrep=nr,
                                   nverb=nv))
   stopifnot(X2.straussm$n == 80)

   # Conditioning upon the number of points of type 1 being 60
   # and the number of points of type 2 being 20:
   X3.straussm <- rmh(model=mod08,start=list(n.start=c(60,20)),
                      control=list(fixall=TRUE,p=1,ptypes=c(0.75,0.25),
                                   nrep=nr,nverb=nv))
   stopifnot(all(table(X3.straussm$marks) == c(60,20)))

   # Multitype Strauss hardcore:
   rhc  <- matrix(c(9.1,5.0,5.0,2.5),2,2)
   mod09 <- list(cif="straushm",par=list(beta=beta,gamma=gmma,
                iradii=r,hradii=rhc),w=c(0,250,0,250))
   X.straushm <- rmh(model=mod09,start=list(n.start=80),
                     control=list(ptypes=c(0.75,0.25),nrep=nr,nverb=nv))

   # Multitype Strauss hardcore with trends for each type:
   beta  <- c(0.27,0.08)
   tr3   <- function(x,y){x <- x/250; y <- y/250;
   			   exp((6*x + 5*y - 18*x^2 + 12*x*y - 9*y^2)/6)
                         }
                         # log quadratic trend
   tr4   <- function(x,y){x <- x/250; y <- y/250;
                         exp(-0.6*x+0.5*y)}
                        # log linear trend
   mod10 <- list(cif="straushm",par=list(beta=beta,gamma=gmma,
                 iradii=r,hradii=rhc),w=c(0,250,0,250),
                 trend=list(tr3,tr4))
   X1.straushm.trend <- rmh(model=mod10,start=list(n.start=350),
                            control=list(ptypes=c(0.75,0.25),
                            nrep=nr,nverb=nv))
   
   # Multitype Strauss hardcore with trends for each type, given as images:
   bigwin <- square(250)
   i1 <- as.im(tr3, bigwin)
   i2 <- as.im(tr4, bigwin)
   mod11 <- list(cif="straushm",par=list(beta=beta,gamma=gmma,
                 iradii=r,hradii=rhc),w=bigwin,
                 trend=list(i1,i2))
   X2.straushm.trend <- rmh(model=mod11,start=list(n.start=350),
                            control=list(ptypes=c(0.75,0.25),expand=1,
                            nrep=nr,nverb=nv))


#######################################################################
############  checks on distribution of output  #######################
#######################################################################

checkp <- function(p, context, testname, failmessage, pcrit=0.01) {
  if(missing(failmessage))
    failmessage <- paste("output failed", testname)
  if(p < pcrit)
    warning(paste(context, ",",  failmessage), call.=FALSE)
  cat(paste("\n", context, ",", testname, "has p-value", signif(p,4), "\n"))
}

# Multitype Strauss code; output is multitype Poisson

beta  <- 100 * c(1,1)
ri    <- matrix(0.07, 2, 2)
gmma  <- matrix(1, 2, 2)  # no interaction
tr1   <- function(x,y){ rep(1, length(x)) }
tr2   <- function(x,y){ rep(2, length(x)) }
mod <- rmhmodel(cif="straussm",
                  par=list(beta=beta,gamma=gmma,radii=ri),
                  w=owin(),
                  trend=list(tr1,tr2))

X <- rmh(mod, start=list(n.start=0), control=list(nrep=1e6))

# The model is Poisson with intensity 100 for type 1 and 200 for type 2.
# Total number of points is Poisson (300)
# Marks are i.i.d. with P(type 1) = 1/3, P(type 2) = 2/3.

# Test whether the total intensity looks right
#
p <- ppois(X$n, 300)
p.val <- 2 * min(p, 1-p)
checkp(p.val, 
       "In multitype Poisson simulation",
       "test whether total number of points has required mean value")

# Test whether the mark distribution looks right
ta <- table(X$marks)
cat("Frequencies of marks:")
print(ta)
checkp(chisq.test(ta, p = c(1,2)/3)$p.value,
       "In multitype Poisson simulation",
       "chi-squared goodness-of-fit test for mark distribution (1/3, 2/3)")

#####
####  multitype Strauss code; fixall=TRUE;
####  output is multinomial process with nonuniform locations
####

the.context <- "In nonuniform multinomial simulation"

beta  <- 100 * c(1,1)
ri    <- matrix(0.07, 2, 2)
gmma  <- matrix(1, 2, 2)  # no interaction
tr1   <- function(x,y){ ifelse(x < 0.5, 0, 2) } 
tr2   <- function(x,y){ ifelse(y < 0.5, 1, 3) }
# cdf of these distributions
Fx1 <- function(x) { ifelse(x < 0.5, 0, ifelse(x < 1, 2 * x - 1, 1)) }
Fy2 <- function(y) { ifelse(y < 0, 0,
                           ifelse(y < 0.5, y/2,
                                  ifelse(y < 1, (1/2 + 3 * (y-1/2))/2, 1))) }
                                                               

mod <- rmhmodel(cif="straussm",
                  par=list(beta=beta,gamma=gmma,radii=ri),
                  w=owin(),
                  trend=list(tr1,tr2))

X <- rmh(mod, start=list(n.start=c(50,50)),
           control=list(nrep=1e6, expand=1, p=1, fixall=TRUE))

# The model is Poisson 
# Mean number of type 1 points = 100
# Mean number of type 2 points = 200
# Total intensity = 300
# Marks are i.i.d. with P(type 1) = 1/3, P(type 2) = 2/3

# Test whether the coordinates look OK
Y <- split(X)
X1 <- Y[[names(Y)[1]]]
X2 <- Y[[names(Y)[2]]]
checkp(ks.test(X1$y, "punif")$p.value,
       the.context,
       "Kolmogorov-Smirnov test of uniformity of y coordinates of type 1 points")
if(any(X1$x < 0.5)) {
  warning(paste(the.context, ",", 
                "x-coordinates of type 1 points are IMPOSSIBLE"), call.=FALSE)
} else {
  checkp(ks.test(Fx1(X1$x), "punif")$p.value,
       the.context,
       "Kolmogorov-Smirnov test of uniformity of transformed x coordinates of type 1 points")
}
checkp(ks.test(X2$x, "punif")$p.value,
       the.context,
     "Kolmogorov-Smirnov test of uniformity of x coordinates of type 2 points")
checkp(ks.test(Fy2(X2$y), "punif")$p.value,
       the.context,
       "Kolmogorov-Smirnov test of uniformity of transformed y coordinates of type 2 points")

})
#
# tests/rmhTrend.R
#
#  Problems with trend images (rmhmodel.ppm or rmhEngine)
#

require(spatstat)
local({
  set.seed(42)

  # Bug folder 37 of 8 feb 2011
  # rmhmodel.ppm -> predict.ppm
  # + rmhResolveTypes -> is.subset.owin

  data(demopat)
  Z <- rescale(demopat, 7000)
  X <- unmark(Z)
  X1 <- split(Z)[[1]]
  Int  <- density(X,dimyx=200)
  Lint <- eval.im(log(npoints(X1)*Int/npoints(X)))
  M    <- as.owin(Int)
  MR   <- intersect.owin(M,scalardilate(M,0.5,origin="midpoint"))
  X1 <- X1[MR]
  Fut  <- ppm(X1,~offset(Lint),covariates=list(Lint=Lint),
              inter=BadGey(r=c(0.03,0.05),sat=3))
  Y   <- rmh(Fut,control=list(expand=M,nrep=1e3), verbose=FALSE)

})
# strange boundary cases

require(spatstat)

local({
   if(!exists("nv"))
     nv <- 0
   if(!exists("nr"))
     nr   <- 5e3

   # Poisson process
   cat("Poisson\n")
   modP <- list(cif="poisson",par=list(beta=10), w = square(3))
   XP <- rmh(model = modP,
             start = list(n.start=25),
             control=list(nrep=nr,nverb=nv))

   # Poisson process case of Strauss
   cat("\nPoisson case of Strauss\n")
   modPS <- list(cif="strauss",par=list(beta=10,gamma=1,r=0.7), w = square(3))
   XPS <- rmh(model=modPS,
              start=list(n.start=25),
              control=list(nrep=nr,nverb=nv))
   
   # Strauss with zero intensity
   cat("\nStrauss with zero intensity\n")
   mod0S <- list(cif="strauss",par=list(beta=0,gamma=0.6,r=0.7), w = square(3))
   X0S   <- rmh(model=mod0S,start=list(n.start=80),
                     control=list(nrep=nr,nverb=nv))
   stopifnot(X0S$n == 0)

   # Poisson with zero intensity
   cat("\nPoisson with zero intensity\n")
   mod0P <- list(cif="poisson",par=list(beta=0), w = square(3))
   X0P <- rmh(model = mod0P,
             start = list(n.start=25),
             control=list(nrep=nr,nverb=nv))


   # Poisson conditioned on zero points
   cat("\nPoisson conditioned on zero points\n")
   modp <- list(cif="poisson",
                 par=list(beta=2), w = square(10))
   Xp <- rmh(modp, start=list(n.start=0), control=list(p=1, nrep=nr))
   stopifnot(Xp$n == 0)

   # Multitype Poisson conditioned on zero points
   cat("\nMultitype Poisson conditioned on zero points\n")
   modp2 <- list(cif="poisson",
                 par=list(beta=2), types=letters[1:3], w = square(10))
   Xp2 <- rmh(modp2, start=list(n.start=0), control=list(p=1, nrep=nr))
   stopifnot(is.marked(Xp2))
   stopifnot(Xp2$n == 0)

   # Multitype Poisson conditioned on zero points of each type
   cat("\nMultitype Poisson conditioned on zero points of each type\n")
   Xp2fix <- rmh(modp2, start=list(n.start=c(0,0,0)),
                 control=list(p=1, fixall=TRUE, nrep=nr))
   stopifnot(is.marked(Xp2fix))
   stopifnot(Xp2fix$n == 0)
    
 })
#
# tests of rmhmodel.ppm
#
require(spatstat)

local({
f <- ppm(cells)
m <- rmhmodel(f)

f <- ppm(cells, ~x)
m <- rmhmodel(f)

f <- ppm(cells, ~1, Strauss(0.1))
m <- rmhmodel(f)

f <- ppm(cells, ~1, StraussHard(r=0.1,hc=0.05))
m <- rmhmodel(f)

f <- ppm(cells, ~1, Hardcore(0.07))
m <- rmhmodel(f)

f <- ppm(cells, ~1, DiggleGratton(0.05,0.1))
m <- rmhmodel(f)

f <- ppm(cells, ~1, Softcore(0.5), correction="isotropic")
m <- rmhmodel(f)

f <- ppm(cells, ~1, Geyer(0.07,2))
m <- rmhmodel(f)

f <- ppm(cells, ~1, BadGey(c(0.07,0.1,0.13),2))
m <- rmhmodel(f)

f <- ppm(cells, ~1, PairPiece(r = c(0.05, 0.1, 0.2)))
m <- rmhmodel(f)

f <- ppm(cells, ~1, AreaInter(r=0.06))
m <- rmhmodel(f)

# multitype

r <- matrix(0.07, 2, 2)
f <- ppm(amacrine, ~1, MultiStrauss(c("off","on"),r))
m <- rmhmodel(f)

h <- matrix(min(nndist(amacrine))/2, 2, 2)
f <- ppm(amacrine, ~1, MultiStraussHard(c("off","on"),r, h))
m <- rmhmodel(f)

diag(r) <- NA
diag(h) <- NA
f <- ppm(amacrine, ~1, MultiStrauss(c("off","on"),r))
m <- rmhmodel(f)

f <- ppm(amacrine, ~1, MultiStraussHard(c("off","on"),r, h))
m <- rmhmodel(f)

# multitype data, interaction not dependent on type

f <- ppm(amacrine, ~marks, Strauss(0.05))
m <- rmhmodel(f)

# trends

f <- ppm(cells, ~x, Strauss(0.1))
m <- rmhmodel(f)

f <- ppm(cells, ~y, StraussHard(r=0.1,hc=0.05))
m <- rmhmodel(f)

f <- ppm(cells, ~x+y, Hardcore(0.07))
m <- rmhmodel(f)

f <- ppm(cells, ~polynom(x,y,2), Softcore(0.5), correction="isotropic")
m <- rmhmodel(f)

# covariates

Z <- as.im(function(x,y){ x^2+y^2 }, as.owin(cells))
f <- ppm(cells, ~z, covariates=list(z=Z))
m <- rmhmodel(f)
m <- rmhmodel(f, control=list(p=1))

Zim <- as.im(Z, as.owin(cells))
f <- ppm(cells, ~z, covariates=list(z=Zim))
m <- rmhmodel(f)

Z <- as.im(function(x,y){ x^2+y }, as.owin(amacrine))
f <- ppm(amacrine, ~z + marks, covariates=list(z=Z))
m <- rmhmodel(f)
m <- rmhmodel(f, control=list(p=1))
m <- rmhmodel(f, control=list(p=1,fixall=TRUE))

Zim <- as.im(Z, as.owin(amacrine))
f <- ppm(amacrine, ~z + marks, covariates=list(z=Zim))
m <- rmhmodel(f)

})
#
#    tests/rmhmodelHybrids.R
#
#  Test that rmhmodel.ppm and rmhmodel.default
#  work on Hybrid interaction models
#
#   $Revision: 1.3 $  $Date: 2013/01/29 02:12:13 $
#

require(spatstat)

local({
  # ......... rmhmodel.ppm .......................

  fit1 <- ppm(redwood, ~1,
              Hybrid(A=Strauss(0.02), B=Geyer(0.1, 2), C=Geyer(0.15, 1)))
  m1 <- rmhmodel(fit1)
  m1
  reach(m1)

  # Test of handling 'IsOffset' 
  data(cells)
  fit2 <- ppm(cells, ~1, Hybrid(H=Hardcore(0.05), G=Geyer(0.15, 2)))
  rmhmodel(fit2)

  # Test of handling Poisson components
  fit3 <- ppm(cells, ~1, Hybrid(P=Poisson(), S=Strauss(0.05)))
  X3 <- rmh(fit3, control=list(nrep=1e3,expand=1), verbose=FALSE)

  # ............ rmhmodel.default ............................

   modH <- list(cif=c("strauss","geyer"),
                par=list(list(beta=50,gamma=0.5, r=0.1),
                         list(beta=1, gamma=0.7, r=0.2, sat=2)),
                w = square(1))
   rmodH <- rmhmodel(modH)
   rmodH
   reach(rmodH)

  # test handling of Poisson components

   modHP <- list(cif=c("poisson","strauss"),
                 par=list(list(beta=5),
                          list(beta=10,gamma=0.5, r=0.1)),
                 w = square(1))
   rmodHP <- rmhmodel(modHP)
   rmodHP
   reach(rmodHP)

   modPP <- list(cif=c("poisson","poisson"),
                 par=list(list(beta=5),
                          list(beta=10)),
                 w = square(1))
   rmodPP <- rmhmodel(modPP)
   rmodPP
   reach(rmodPP)
  
})


#
#  tests/rmh.ppm.R
#
#  $Revision: 1.1 $ $Date: 2012/10/14 07:24:21 $
#
#  Examples removed from rmh.ppm.Rd
#  stripped down to minimal tests of validity
#

require(spatstat)
local({
   op <- spatstat.options()
   spatstat.options(rmh.nrep=10, npixel=10, ndummy.min=10)
   spatstat.options(project.fast=TRUE)
   Nrep <- 10

   X <- swedishpines
   # Poisson process
   fit <- ppm(X, ~1, Poisson())
   Xsim <- rmh(fit)
   # Strauss process   
   fit <- ppm(X, ~1, Strauss(r=7))
   Xsim <- rmh(fit)

   # Strauss process simulated on a larger window
   # then clipped to original window
   Xsim <- rmh(fit, control=list(nrep=Nrep, expand=1.1, periodic=TRUE))

   # Strauss - hard core process
#   fit <- ppm(X, ~1, StraussHard(r=7,hc=2))
#   Xsim <- rmh(fit, start=list(n.start=X$n))

   # Geyer saturation process
#   fit <- ppm(X, ~1, Geyer(r=7,sat=2))
#   Xsim <- rmh(fit, start=list(n.start=X$n))

   # Area-interaction process
     fit <- ppm(X, ~1, AreaInter(r=7))
     Xsim <- rmh(fit, start=list(n.start=X$n))
  
     # soft core interaction process
#     X <- quadscheme(X, nd=50)
#     fit <- ppm(X, ~1, Softcore(kappa=0.1), correction="isotropic")
#     Xsim <- rmh(fit, start=list(n.start=X$n))

     # Diggle-Gratton pairwise interaction model
#     fit <- ppm(cells, ~1, DiggleGratton(0.05, 0.1))
#     Xsim <- rmh(fit, start=list(n.start=cells$n))
#     plot(Xsim, main="simulation from fitted Diggle-Gratton model")
   
   X <- rSSI(0.05, 100)

   # piecewise-constant pairwise interaction function
   fit <- ppm(X, ~1, PairPiece(seq(0.02, 0.1, by=0.01)))
   Xsim <- rmh(fit)

   # marked point pattern
   Y <- amacrine

   # marked Poisson models
   fit <- ppm(Y)
   Ysim <- rmh(fit)

   fit <- ppm(Y,~marks)
   Ysim <- rmh(fit)

   fit <- ppm(Y,~x)
   Ysim <- rmh(fit)
#   fit <- ppm(Y,~polynom(x,2))
#   Ysim <- rmh(fit)

   fit <- ppm(Y,~marks+x)
   Ysim <- rmh(fit)
#   fit <- ppm(Y,~marks+polynom(x,2))
#   Ysim <- rmh(fit)

   # multitype Strauss models
   MS <- MultiStrauss(types = levels(Y$marks),
                      radii=matrix(0.07, ncol=2, nrow=2))

#   fit <- ppm(Y,~marks*polynom(x,2), MS)
    fit <- ppm(Y,~marks*x, MS)
   Ysim <- rmh(fit)

   spatstat.options(op)
 })
# fvproblems.R

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

  levels(marks(amacrine)) <- c("Nastricreechia krorluppia", "Homo habilis")
  plot(Kcross(amacrine))
  plot(alltypes(amacrine, "K"))
  plot(alltypes(amacrine, "J"))
  plot(alltypes(amacrine, pcfcross))
})
# test of case where mark levels contain illegal characters

require(spatstat)
local({
  hyphenated <- c("a", "not-a")
  spaced <- c("U", "non U")
  suffixed <- c("a+", "a*")
  charred <- c("+", "*")

  irad <- matrix(0.1, 2,2)
  hrad <- matrix(0.005, 2, 2)

  tryit <- function(types, X, irad, hrad) { 
    levels(marks(X)) <- types
    fit <- ppm(X, ~marks + polynom(x,y,2),
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
})
#
# test cases where there are no (rows or columns of) marks
#

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
# checks validity of fast C implementation of Geyer interaction
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
})

require(spatstat)
local({
  co <- as.ppp(corners(letterR), letterR, check=FALSE)
  co[letterR]
})

#  tests/segments.R
#  $Revision: 1.7 $  $Date: 2011/12/05 07:29:16 $

require(spatstat)

local({
# pointed out by Jeff Laake
W <- owin()
X <- psp(x0=.25,x1=.25,y0=0,y1=1,window=W)
X[W]

# migrated from 'lpp'

X <- psp(runif(10),runif(10),runif(10),runif(10), window=owin())
Z <- as.mask.psp(X)
Z <- pixellate(X)

# test of distppll pointed out by Ang Qi Wei

p <- matrix(c(1.5, 0), 1, 2)
l <- matrix(c(0,0,1,0,1,0,2,0), 2, 4, byrow=T)
a <- distppll(p, l, mintype=2, method="interpreted")
b <- distppll(p, l, mintype=2, method="Fortran")
d <- distppll(p, l, mintype=2, method="C")
if(a$min.which != b$min.which)
  stop("conflict between Fortran and interpreted code in distppll")
if(a$min.which != d$min.which)
  stop("conflict between C and interpreted code in distppll")

# tests of pixellate.psp -> seg2pixL

ns <- 50
out <- numeric(ns)
for(i in 1:ns) {
  X <- psp(runif(1), runif(1), runif(1), runif(1), window=owin())
  len <- lengths.psp(X)
  dlen <- sum(pixellate(X)$v)
  out[i] <- if(len > 1e-7) dlen/len else 1
}
if(diff(range(out)) > 0.01) stop(paste(
       "pixellate.psp test 1: relative error [",
       paste(diff(range(out)), collapse=", "),
       "]"))

# Michael Sumner's test examples

set.seed(33)
n <- 2001
co <- cbind(runif(n), runif(n))
ow <- owin()
X <- psp(co[-n,1], co[-n,2], co[-1,1], co[-1,2], window=ow)
s1 <- sum(pixellate(X))
s2 <- sum(lengths.psp(X))
if(abs(s1 - s2)/s2 > 0.01) {
  stop(paste("pixellate.psp test 2:",
             "sum(pixellate(X)) = ", s1,
             "!=", s2, "= sum(lengths.psp(X))"))
}

wts <- 1/(lengths.psp(X) * X$n)
s1 <- sum(pixellate(X, weights=wts))
if(abs(s1-1) > 0.01) {
  stop(paste("pixellate.psp test 3:",
             "sum(pixellate(X, weights))=", s1,
             " (should be 1)"))
}

X <- psp(0, 0, 0.01, 0.001, window=owin())
s1 <- sum(pixellate(X))
s2 <- sum(lengths.psp(X))
if(abs(s1 - s2)/s2 > 0.01) {
  stop(paste("pixellate.psp test 4:",
             "sum(pixellate(X)) = ", s1,
             "!=", s2, "= sum(lengths.psp(X))"))
}

X <- psp(0, 0, 0.001, 0.001, window=owin())
s1 <- sum(pixellate(X))
s2 <- sum(lengths.psp(X))
if(abs(s1 - s2)/s2 > 0.01) {
  stop(paste("pixellate.psp test 5:",
             "sum(pixellate(X)) = ", s1,
             "!=", s2, "= sum(lengths.psp(X))"))
}

})
# test for step() operation
#
require(spatstat)
local({
  Z <- as.im(function(x,y){ x^3 - y^2 }, nztrees$window)
  fitP <- ppm(nztrees, ~x+y+Z, covariates=list(Z=Z))
  step(fitP)
  fitS <- update(fitP, Strauss(7))
  step(fitS)
  fitM <- ppm(amacrine, ~ marks*(x+y),
              MultiStrauss(types=levels(marks(amacrine)), radii=matrix(0.04, 2, 2)))
  step(fitM)
})


require(spatstat)
source("badwindow.R")
owinpolycheck(W,verbose=FALSE)




#       
#        tests/hobjects.R
#
#   Validity of methods for ppm(... method="ho")
#

require(spatstat)

local({
  fit  <- ppm(cells, ~1,         Strauss(0.1), method="ho", nsim=10)
  fitx <- ppm(cells, ~offset(x), Strauss(0.1), method="ho", nsim=10)

  a  <- AIC(fit)
  ax <- AIC(fitx)

  f  <- fitted(fit)
  fx <- fitted(fitx)

  p  <- predict(fit)
  px <- predict(fitx)
})


# check kstest with strange data
require(spatstat)
local({
  # Marked point patterns with some marks not represented
  AC <- split(ants, un=FALSE)$Cataglyphis
  AM <- split(ants, un=FALSE)$Messor
  DM <- distmap(AM)
  # should produce a warning, rather than a crash:
  kstest(AC, DM)
  # should be OK:
  kstest(unmark(AC), DM)
  # linear networks
  X <- runiflpp(20, simplenet)
  fit <- lppm(X, ~1)
  kstest(fit, "y")
})

#
# tests/kppm.R
#
# $Revision: 1.6 $ $Date: 2012/04/08 03:22:20 $
#
# Test functionality of kppm that depends on RandomFields
#

require(spatstat)
local({

  if(require(RandomFields) && RandomFieldsSafe()) {

    fit0 <- kppm(redwood, ~1, "LGCP")
    simulate(fit0)

    fit <- kppm(redwood, ~x, "LGCP",
                covmodel=list(model="matern", nu=0.3),
                control=list(maxit=5))
    simulate(fit)

# ... and Abdollah's code

    fit <- kppm(redwood, ~x, cluster="Cauchy", statistic="K")
    simulate(fit)
  }
  
})


# temporary test file for localpcfmatrix

require(spatstat)
local({
  a <- localpcfmatrix(redwood)
  a
  plot(a)
  a[, 3:5]
})
# check for various bugs related to factor conversions
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
})


#
#  tests/splitpea.R
#
#  Check behaviour of split.ppp etc
#
#  Thanks to Marcelino de la Cruz
#
#  $Revision: 1.8 $  $Date: 2013/10/06 08:44:28 $
#

require(spatstat)

local({
W <- square(8)
X <- ppp(c(2.98, 4.58, 7.27, 1.61, 7.19),
         c(7.56, 5.29, 5.03, 0.49, 1.65),
         window=W)
Z <- quadrats(W, 4, 4)
Yall <- split(X, Z, drop=FALSE)
Ydrop <- split(X, Z, drop=TRUE)

P <- Yall[[1]]
if(!all(inside.owin(P$x, P$y, P$window)))
  stop("Black hole detected when drop=FALSE")
P <- Ydrop[[1]]
if(!all(inside.owin(P$x, P$y, P$window)))
  stop("Black hole detected when drop=TRUE")

Ydrop[[1]] <- P[1]
split(X, Z, drop=TRUE) <- Ydrop

# test NA handling
Zbad <- quadrats(square(4), 2, 2)
Ybdrop <- split(X, Zbad, drop=TRUE)
Yball  <- split(X, Zbad, drop=FALSE)

# From Marcelino
set.seed(1)
W<- square(10) # the big window
puntos<- rpoispp(0.5, win=W)
data(letterR)
r00 <- letterR
r05 <- shift(letterR,c(0,5))
r50 <- shift(letterR,c(5,0))
r55 <- shift(letterR,c(5,5))
tessr4 <- tess(tiles=list(r00, r05,r50,r55))
puntosr4 <- split(puntos, tessr4, drop=TRUE)
split(puntos, tessr4, drop=TRUE) <- puntosr4

})
#
#  tests/imageops.R
#
#   $Revision: 1.5 $   $Date: 2011/12/05 07:29:16 $
#

require(spatstat)
local({
  A <- as.im(owin())
  B <- as.im(owin(c(1.1, 1.9), c(0,1)))
  Z <- imcov(A, B)
  stopifnot(abs(max(Z) - 0.8) < 0.1)
})




# tests/triplets.R
# test code for triplet interaction
# $Revision: 1.4 $ $Date: 2012/07/12 02:43:32 $
require(spatstat)
local({
  fit <- ppm(redwood, ~1, Triplets(0.1))
  fit
  suffstat(fit)
  # hard core (zero triangles, coefficient is NA)
  fit0 <- ppm(cells, ~1, Triplets(0.05))
  fit0
  suffstat(fit0)
  # bug case (1 triangle in data)
  fit1 <- ppm(cells, ~1, Triplets(0.15))
  fit1
  suffstat(fit1)
})
#
#     tests/project.ppm.R
#
#      $Revision: 1.3 $  $Date: 2012/10/22 03:12:08 $
#
#     Tests of projection mechanism
#

require(spatstat)
local({
  # a very unidentifiable model
  fit <- ppm(cells, ~Z, Strauss(1e-06), covariates=list(Z=0))
  project.ppm(fit)
  # multitype
  fit2 <- ppm(amacrine, ~1, MultiStrauss(types=c("off", "on"),
                                         radii=matrix(1e-06, 2, 2)))
  project.ppm(fit2)
  
  # hybrids
  r0 <- min(nndist(redwood))
  ra <- 1.25 * r0
  rb <- 0.8 * r0
  f <- ppm(redwood, ~1, Hybrid(A=Strauss(ra), B=Geyer(0.1, 2)), project=TRUE)
  f <- ppm(redwood, ~1, Hybrid(A=Strauss(rb), B=Geyer(0.1, 2)), project=TRUE)
  f <- ppm(redwood, ~1, Hybrid(A=Strauss(ra), B=Strauss(0.1)), project=TRUE)
  f <- ppm(redwood, ~1, Hybrid(A=Strauss(rb), B=Strauss(0.1)), project=TRUE)
  f <- ppm(redwood, ~1, Hybrid(A=Hardcore(rb), B=Strauss(0.1)), project=TRUE)
  f <- ppm(redwood, ~1, Hybrid(A=Hardcore(rb), B=Geyer(0.1, 2)), project=TRUE)
  f <- ppm(redwood, ~1, Hybrid(A=Geyer(rb, 1), B=Strauss(0.1)), project=TRUE)

})
#
# tests/hyperframe.R
#
# test "[.hyperframe" etc
#
#  $Revision: 1.2 $  $Date: 2012/01/31 11:04:44 $
#

  lambda <- runif(4, min=50, max=100)
  X <- lapply(as.list(lambda), function(x) { rpoispp(x) })
  h <- hyperframe(lambda=lambda, X=X)
  h$lambda2 <- lambda^2
  h[, "lambda3"] <- lambda^3
  h[, "Y"] <- X
  h[, "X"] <- lapply(X, flipxy)
  h[, c("X", "Y")] <- hyperframe(X=X, Y=X)

# check fast code for Kest
require(spatstat)
local({
  Kb <- Kest(cells, nlarge=0)
  Ku <- Kest(cells, correction="none")
  Kbu <- Kest(cells, correction=c("none", "border"))
})


#
#  tests/vcovppm.R
#
#  Check validity of vcov.ppm algorithms
#
#  Thanks to Ege Rubak
#
#  $Revision: 1.4 $  $Date: 2013/09/20 09:01:34 $
#

require(spatstat)

local({

  set.seed(42)
  X <- rStrauss(200, .5, .05)
  model <- ppm(X, inter = Strauss(.05))

  b  <- vcov(model, generic = TRUE, algorithm = "basic")
  v  <- vcov(model, generic = TRUE, algorithm = "vector")
  vc <- vcov(model, generic = TRUE, algorithm = "vectorclip")
  vn <- vcov(model, generic = FALSE)

  disagree <- function(x, y, tol=1e-7) { max(abs(x-y)) > tol }
  asymmetric <- function(x) { disagree(x, t(x)) }

  if(asymmetric(b))
    stop("Non-symmetric matrix produced by vcov.ppm 'basic' algorithm")
  if(asymmetric(v))
    stop("Non-symmetric matrix produced by vcov.ppm 'vector' algorithm")
  if(asymmetric(vc))
    stop("Non-symmetric matrix produced by vcov.ppm 'vectorclip' algorithm")
  if(asymmetric(vn))
    stop("Non-symmetric matrix produced by vcov.ppm Strauss algorithm")
  
  if(disagree(v, b))
    stop("Disagreement between vcov.ppm algorithms 'vector' and 'basic' ")
  if(disagree(v, vc))
    stop("Disagreement between vcov.ppm algorithms 'vector' and 'vectorclip' ")
  if(disagree(vn, vc))
    stop("Disagreement between vcov.ppm generic and Strauss algorithms")

  # Geyer code
  xx <- c(0.7375956, 0.6851697, 0.6399788, 0.6188382)
  yy <- c(0.5816040, 0.6456319, 0.5150633, 0.6191592)
  Y <- ppp(xx, yy, window=square(1))
  modelY <- ppm(Y, ~1, Geyer(0.1, 1))

  b  <- vcov(modelY, generic = TRUE, algorithm = "basic")
  v  <- vcov(modelY, generic = TRUE, algorithm = "vector")
  vc <- vcov(modelY, generic = TRUE, algorithm = "vectorclip")

  if(asymmetric(b))
    stop("Non-symmetric matrix produced by vcov.ppm 'basic' algorithm for Geyer model")
  if(asymmetric(v))
    stop("Non-symmetric matrix produced by vcov.ppm 'vector' algorithm for Geyer model")
  if(asymmetric(vc))
    stop("Non-symmetric matrix produced by vcov.ppm 'vectorclip' algorithm for Geyer model")
  
  if(disagree(v, b))
    stop("Disagreement between vcov.ppm algorithms 'vector' and 'basic' for Geyer model")
  if(disagree(v, vc))
    stop("Disagreement between vcov.ppm algorithms 'vector' and 'vectorclip' for Geyer model")

  # test of pairwise.family$delta2
  modelZ <- ppm(amacrine, ~1, MultiStrauss(radii=matrix(0.1, 2, 2)))
  b <- vcov(modelZ)
  g <- vcov(modelZ, generic=TRUE)
  if(disagree(b, g))
    stop("Disagreement between vcov.ppm algorithms for MultiStrauss model")
  
})
# tests/windows.R
# Tests of owin geometry code

require(spatstat)
local({
  # Ege Rubak spotted this problem in 1.28-1
  A <- as.owin(ants)
  B <- dilation(A, 140)
  if(!is.subset.owin(A, B))
    stop("is.subset.owin fails in polygonal case")

  # thanks to Tom Rosenbaum
  A <- shift(square(3), origin="midpoint")
  B <- shift(square(1), origin="midpoint")
  AB <- setminus.owin(A, B)
  D <- shift(square(2), origin="midpoint")
  if(is.subset.owin(D,AB))
    stop("is.subset.owin fails for polygons with holes")
})

#
# test addvar options
#

X <-  rpoispp(function(x,y){exp(3+3*x)})
model <- ppm(X, ~y)
addvar(model, "x", crosscheck=TRUE)
addvar(model, "x", bw.input="quad")
w <- square(0.5)
addvar(model, "x", subregion=w)
addvar(model, "x", subregion=w, bw.input="points")
# additional test of parres
X <-  rpoispp(function(x,y){exp(3+x+2*x^2)})
model <- ppm(X, ~x+y)

# options in parres
parres(model, "x")
parres(model, "x", bw.input="quad")
w <- square(0.5)
parres(model, "x", subregion=w)
parres(model, "x", subregion=w, bw.input="quad")

# check whether 'update.ppm' has messed up internals
mod2 <- update(model, ~x)
parres(mod2, "x")
#
# tests/lppstuff.R
#
# Tests for lpp code
#
#  $Revision: 1.2 $  $Date: 2013/01/23 08:20:41 $


require(spatstat)

local({
  # check 'normalise' option in linearKinhom
  X <- rpoislpp(5, simplenet)
  fit <- lppm(X, ~x)
  K <- linearKinhom(X, lambda=fit, normalise=FALSE)
  plot(K)
  g <- linearpcfinhom(X, lambda=fit, normalise=FALSE)
  plot(g)
  K <- linearKinhom(X, lambda=fit, normalise=TRUE)
  plot(K)
  g <- linearpcfinhom(X, lambda=fit, normalise=TRUE)
  plot(g)
  # check empty patterns OK
  X <- runiflpp(0, simplenet)
  print(X)
})

#
#  tests/density.R
#
#  Test behaviour of density methods and inhomogeneous summary functions
#
#  $Revision: 1.1 $  $Date: 2013/02/26 09:13:52 $
#

require(spatstat)

local({

  lam <- density(redwood)
  K <- Kinhom(redwood, lam)
  
  lamX <- density(redwood, at="points")
  KX <- Kinhom(redwood, lamX)
  
})
#
# tests/slrm.R
#
# $Revision: 1.1 $ $Date: 2013/04/19 10:14:52 $
#
# Test slrm fitting and prediction when there are NA's
#

require(spatstat)
local({
  X <- copper$SouthPoints
  W <- owin(poly=list(x=c(0,35,35,1),y=c(1,1,150,150)))
  Y <- X[W]
  fit <- slrm(Y ~ x+y)
  pred <- predict(fit)
})


# tests/linalgeb.R
# checks validity of linear algebra code
#  $Revision: 1.2 $ $Date: 2013/04/18 09:14:37 $
require(spatstat)
local({
  p <- 3
  n <- 4
  x <- array(as.numeric(1:(p * n * n)), dim=c(p, n, n))
  w <- matrix(1:(n*n), n, n)
  y <- matrix(numeric(p * p), p, p)
  for(i in 1:n)
    for(j in (1:n)[-i])
      y <- y + w[i,j] * outer(x[,i,j], x[,j,i])
  z <- sumsymouter(x, w)
  if(!identical(y,z))
    stop("sumsymouter gives incorrect result")
})
#
#  tests/undoc.R
#
#   $Revision: 1.1 $   $Date: 2013/07/25 10:26:09 $
#
#  Test undocumented hacks, etc

require(spatstat)
local({
  # pixellate.ppp accepts a data frame of weights
  pixellate(cells, weights=data.frame(a=1:42, b=42:1))
})




#
# tests/mppm.R
#
# Basic tests of mppm
#
# $Revision: 1.1 $ $Date: 2013/11/10 08:59:08 $
# 

require(spatstat)

local({
# test interaction formulae and subfits
fit1 <- mppm(Points ~ group, simba, hyperframe(po=Poisson(), str=Strauss(0.1)),
            iformula=~ifelse(group=="control", po, str))
fit2 <- mppm(Points ~ group, simba, hyperframe(po=Poisson(), str=Strauss(0.1)),
            iformula=~id * str)
fit3 <- mppm(Points ~ group, simba, hyperframe(po=Poisson(), pie=PairPiece(c(0.05,0.1))), iformula=~I((group=="control") * po) + I((group=="treatment") * pie))
fit1
fit2
fit3
subfits(fit1)
subfits(fit2)
subfits(fit3)

# test handling of offsets and zero cif values in mppm

 data(waterstriders)
 H <- hyperframe(Y = waterstriders)
 mppm(Y ~ 1,  data=H, Hardcore(1.5))
 mppm(Y ~ 1,  data=H, StraussHard(7, 1.5))
})
# tests/ppx.R
#
# Test operations for ppx objects
#
#  $Revision: 1.1 $ $Date: 2013/11/19 03:36:27 $
#

require(spatstat)

local({
  df <- data.frame(x=c(1,2,2,1), y=c(1,2,3,1), z=c(2,3,4,2))
  X <- ppx(data=df, coord.type=rep("s", 3), domain=box3())
  unique(X)
  duplicated(X)
  multiplicity(X)
})
