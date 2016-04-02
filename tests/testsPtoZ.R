## 
##    tests/percy.R
##
## Tests of Percus-Yevick approximations
##
##    $Revision: 1.2 $ $Date: 2015/12/29 08:54:49 $

require(spatstat)
local({
  fit <- ppm(swedishpines ~1, DiggleGatesStibbard(6))
  K <- Kmodel(fit)
})

##
## tests/pixelgripes.R
##     Problems related to pixellation of windows
##
## $Revision: 1.3 $ $Date: 2015/12/29 08:54:49 $

require(spatstat)
local({
  ## From Philipp Hunziker: bug in rNeymanScott (etc)
  ## Create an irregular window
  PM <- matrix(c(1,0,0.5,1,0,0), 3, 2, byrow=TRUE)
  P <- owin(poly=PM)
  ## Generate Matern points
  X <- rMatClust(50, 0.05, 5, win=P)
  ## Some distance function as a covariate
  distorigin <- function(x, y) { sqrt(x^2 + y^2) }
  ## No covariates: works fine
  fit0 <- kppm(X ~ 1, clusters="MatClust")
  Y0 <- simulate(fit0, retry=0)
  ## Covariates: Simulation fails
  fit1 <- kppm(X ~ distorigin, clusters="MatClust")
  Y1 <- simulate(fit1, retry=0)
})
## 
## tests/polygons.R
##
##  $Revision: 1.2 $ $Date: 2015/12/29 08:54:49 $
##
require(spatstat)
local({
  co <- as.ppp(corners(letterR), letterR, check=FALSE)
  co[letterR]
})

# 
#   tests/ppmBadData.R
#
# $Revision: 1.5 $ $Date: 2015/12/29 08:54:49 $

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
  fY <- ppm(Y ~cv + I(cv^2), data=list(cv=M), correction="translate")
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
  ljtrmod <- ppm(ljtr, trend= ~ Z, interaction=NULL, data=list(Z=rain))
  ss <- summary(ljtrmod)
})

local({
  # From Ege
  # Degenerate but non-null argument 'covariates'
  xx <- list()
  names(xx) <- character(0)
  fit <- ppm(cells ~x, covariates = xx)
  st <- summary(fit) 
})
#
#   tests/ppmgam.R
#
#   Test ppm with use.gam=TRUE
#
#   $Revision: 1.3 $  $Date: 2015/09/01 02:01:33 $
#

require(spatstat)
local({
  fit <- ppm(nztrees ~s(x,y), use.gam=TRUE)
  mm <- model.matrix(fit)
  mf <- model.frame(fit)
  v <- vcov(fit)
  prd <- predict(fit)
})

#
#  tests/ppmlogi.R
#
# Tests of ppm(method='logi')
#
# $Revision: 1.3 $  Date$
#

require(spatstat)
local({
  fit <- ppm(cells ~x, method="logi")
  f <- fitted(fit)
  p <- predict(fit)
  fitS <- ppm(cells ~x, Strauss(0.08), method="logi")
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
#
#   tests/ppmmarkorder.R
#
# $Revision: 1.3 $  $Date: 2015/12/29 08:54:49 $
#
# Test that predict.ppm, plot.ppm and plot.fitin
# tolerate marks with levels that are not in alpha order
#
require(spatstat)
local({
  X <- amacrine
  levels(marks(X)) <- c("ZZZ", "AAA")
  fit <- ppm(X ~marks, MultiStrauss(c("ZZZ","AAA"), matrix(0.06, 2, 2)))
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
#   $Revision: 1.5 $  $Date: 2015/12/29 08:54:49 $
#

require(spatstat)
local({
  ##   (1) Scoping problem that can arise when ppm splits the data

  fit <- ppm(bei ~elev, data=bei.extra)
  mm <- model.matrix(fit)

  ##   (2) Fast update mechanism

  fit1 <- ppm(cells ~x+y, Strauss(0.07))
  fit2 <- update(fit1, ~y)
  fit3 <- update(fit2, ~x)

  ## (3) New formula-based syntax
  attach(bei.extra)
  slfit <- ppm(bei ~ grad)
  sl2fit <- update(slfit, ~grad + I(grad^2))
  slfitup <- update(slfit, use.internal=TRUE)
  sl2fitup <- update(sl2fit, use.internal=TRUE)
})

#
#   tests/ppmtricks.R
#
#   Test backdoor exits and hidden options in ppm
#        and summary.ppm, print.summary.ppm
#
#   $Revision: 1.5 $  $Date: 2015/12/29 08:54:49 $
#
require(spatstat)
local({

  ## (1) skip.border
  
  fit <- ppm(cells, ~1, Strauss(0.1), skip.border=TRUE)

  ## (2) subset arguments of different kinds
  fut <- ppm(cells ~ x, subset=(x > 0.5))
  fot <- ppm(cells ~ x, subset=(x > 0.5), method="logi")
  W <- owin(c(0.4, 0.8), c(0.2, 0.7))
  fut <- ppm(cells ~ x, subset=W)
  fot <- ppm(cells ~ x, subset=W, method="logi")

  ## (3) profilepl -> ppm
  ##     uses 'skip.border' and 'precomputed'
  ##     also tests scoping for covariates
  splants <- split(ants)
  mess    <- splants[["Messor"]]
  cats    <- splants[["Cataglyphis"]]
  ss      <- data.frame(r=seq(60,120,by=20),hc=29/6)
  dM      <- distmap(mess,dimyx=256)
  mungf    <- profilepl(ss, StraussHard, cats ~ dM)
  mungp   <- profilepl(ss, StraussHard, trend=~dM, Q=cats)

  ## (4) splitting large quadschemes into blocks
  op <- spatstat.options(maxmatrix=5000)
  pr <- predict(ppm(cells ~ x, AreaInter(0.05)))

  ## (5) shortcuts in summary.ppm
  ## and corresponding behaviour of print.summary.ppm
  print(summary(fit, quick=TRUE))
  print(summary(fit, quick="entries"))
  print(summary(fit, quick="no prediction"))
  print(summary(fit, quick="no variances"))
  
  spatstat.options(op)
})

#
# tests/ppx.R
#
# Test operations for ppx objects
#
#  $Revision: 1.2 $ $Date: 2015/12/29 08:54:49 $
#

require(spatstat)

local({
  df <- data.frame(x=c(1,2,2,1), y=c(1,2,3,1), z=c(2,3,4,2))
  X <- ppx(data=df, coord.type=rep("s", 3), domain=box3())
  unique(X)
  duplicated(X)
  multiplicity(X)
})
#
# tests/prediction.R
#
# Things that might go wrong with predict()
#
#  $Revision: 1.4 $ $Date: 2016/03/04 03:14:40 $
#

require(spatstat)

local({
  # test of 'covfunargs'
  f <- function(x,y,a){ y - a }
  fit <- ppm(cells ~x + f, covariates=list(f=f), covfunargs=list(a=1/2))
  p <- predict(fit)

  # prediction involving 0 * NA
  qc <- quadscheme(cells, nd=10)
  r <- minnndist(as.ppp(qc))/10
  fit <- ppm(qc ~ 1, Strauss(r)) # model has NA for interaction coefficient
  p1 <- predict(fit)
  p2 <- predict(fit, type="cif", ngrid=10)
  stopifnot(all(is.finite(as.matrix(p1))))
  stopifnot(all(is.finite(as.matrix(p2))))

  # test of 'new.coef' mechanism
  fut <- ppm(cells ~ x, Strauss(0.15), rbord=0)
  p0 <- predict(fut, type="cif")
  pe <- predict(fut, type="cif", new.coef=coef(fut))
  pn <- predict(fut, type="cif", new.coef=unname(coef(fut)))
  if(max(abs(pe-p0)) > 0.01)
    stop("new.coef mechanism is broken!")
  if(max(abs(pn-p0)) > 0.01)
    stop("new.coef mechanism gives wrong answer, for unnamed vectors")

  # tests of relrisk.ppm
  fut <- ppm(amacrine ~ x * marks)
  a <- relrisk(fut, control=2, relative=TRUE)
  a <- relrisk(fut, se=TRUE)
  a <- relrisk(fut, relative=TRUE, se=TRUE)
  fut <- ppm(sporophores ~ marks + x)
  a <- relrisk(fut, control=2, relative=TRUE)
  a <- relrisk(fut, se=TRUE)
  a <- relrisk(fut, relative=TRUE, se=TRUE)
  
})

#
#     tests/project.ppm.R
#
#      $Revision: 1.6 $  $Date: 2015/08/27 08:19:03 $
#
#     Tests of projection mechanism
#

require(spatstat)
local({
  chk <- function(m) {
    if(!valid.ppm(m)) stop("Projected model was still not valid")
    return(invisible(NULL))
  }
  # a very unidentifiable model
  fit <- ppm(cells ~Z, Strauss(1e-06), covariates=list(Z=0))
  chk(emend(fit))
  # multitype
  r <- matrix(1e-06, 2, 2)
  fit2 <- ppm(amacrine ~1, MultiStrauss(types=c("off", "on"), radii=r))
  chk(emend(fit2))
  # complicated multitype 
  fit3 <- ppm(amacrine ~1, MultiStraussHard(types=c("off", "on"),
                                            iradii=r, hradii=r/5))
  chk(emend(fit3))
  
  # hierarchical
  ra <- r
  r[2,1] <- NA
  fit4 <- ppm(amacrine ~1, HierStrauss(types=c("off", "on"), radii=r))
  chk(emend(fit4))
  # complicated hierarchical
  fit5 <- ppm(amacrine ~1, HierStraussHard(types=c("off", "on"),
                                            iradii=r, hradii=r/5))
  chk(emend(fit5))
  
  # hybrids
  r0 <- min(nndist(redwood))
  ra <- 1.25 * r0
  rb <- 0.8 * r0
  f1 <- ppm(redwood ~1, Hybrid(A=Strauss(ra), B=Geyer(0.1, 2)), project=TRUE)
  chk(f1)
  f2 <- ppm(redwood ~1, Hybrid(A=Strauss(rb), B=Geyer(0.1, 2)), project=TRUE)
  chk(f2)
  f3 <- ppm(redwood ~1, Hybrid(A=Strauss(ra), B=Strauss(0.1)), project=TRUE)
  chk(f3)
  f4 <- ppm(redwood ~1, Hybrid(A=Strauss(rb), B=Strauss(0.1)), project=TRUE)
  chk(f4)
  f5 <- ppm(redwood ~1, Hybrid(A=Hardcore(rb), B=Strauss(0.1)), project=TRUE)
  chk(f5)
  f6 <- ppm(redwood ~1, Hybrid(A=Hardcore(rb), B=Geyer(0.1, 2)), project=TRUE)
  chk(f6)
  f7 <- ppm(redwood ~1, Hybrid(A=Geyer(rb, 1), B=Strauss(0.1)), project=TRUE)
  chk(f7)

})
#'  tests/randoms.R
#'   Further tests of random generation code
#'  $Revision: 1.1 $ $Date: 2016/04/02 04:03:46 $

require(spatstat)
local({
  A <- runifrect(6, nsim=2)
  A <- runifdisc(6, nsim=2)
  A <- runifpoispp(5, nsim=2)
  A <- runifpoispp(0, nsim=2)
  Z <- as.im(function(x,y) 10*x, square(1))
  A <- rpoint(n=6, f=Z, fmax=10, nsim=2)
  A <- rSSI(0.05, 6, nsim=2)
  A <- rstrat(nx=4, nsim=2)
  A <- rsyst(nx=4, nsim=2)
  A <- rthin(cells, P=0.5, nsim=2)
  A <- rjitter(cells, nsim=2, retry=FALSE)
})

##
## tests/rhohat.R
##
## Test all combinations of options for rhohatCalc
##
## $Revision: 1.2 $ $Date: 2015/12/29 08:54:49 $

local({
  require(spatstat)
  X <-  rpoispp(function(x,y){exp(3+3*x)})
  ## done in example(rhohat):
  ## rhoA <- rhohat(X, "x")
  ## rhoB <- rhohat(X, "x", method="reweight")
  ## rhoC <- rhohat(X, "x", method="transform")
  fit <- ppm(X, ~x)
  rhofitA <- rhohat(fit, "x")
  rhofitB <- rhohat(fit, "x", method="reweight")
  rhofitC <- rhohat(fit, "x", method="transform")

  ## Baseline
  lam <- predict(fit)
  rhoAb <- rhohat(X, "x", baseline=lam)
  rhoBb <- rhohat(X, "x", method="reweight", baseline=lam)
  rhoCb <- rhohat(X, "x", method="transform", baseline=lam)

  ## Horvitz-Thompson
  rhoAH <- rhohat(X, "x", horvitz=TRUE) 
  rhoBH <- rhohat(X, "x", method="reweight", horvitz=TRUE)
  rhoCH <- rhohat(X, "x", method="transform", horvitz=TRUE)
  rhofitAH <- rhohat(fit, "x", horvitz=TRUE)
  rhofitBH <- rhohat(fit, "x", method="reweight", horvitz=TRUE)
  rhofitCH <- rhohat(fit, "x", method="transform", horvitz=TRUE)
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
##
##   tests/rmhBasic.R
##
##   $Revision: 1.11 $  $Date: 2015/12/29 08:54:49 $
#
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
##
##     tests/rmhErrors.R
##
##     $Revision: 1.5 $ $Date: 2015/12/29 08:54:49 $
##
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
#  tests/rmhMulti.R
#
#  tests of rmh, running multitype point processes
#
#   $Revision: 1.6 $  $Date: 2015/12/29 08:54:49 $

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
  stop(paste(the.context, ",", 
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
  Fut  <- ppm(X1~offset(Lint),covariates=list(Lint=Lint),
              inter=BadGey(r=c(0.03,0.05),sat=3))
  Y   <- rmh(Fut,control=list(expand=M,nrep=1e3), verbose=FALSE)

})
#
#   tests/rmhWeird.R
#
#   $Revision: 1.3 $  $Date: 2015/12/29 08:54:49 $
#
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
#      tests/rmhmodel.ppm.R
#
#    $Revision: 1.8 $  $Date: 2015/12/29 08:54:49 $
#
# Case-by-case tests of rmhmodel.ppm
#
require(spatstat)

local({
f <- ppm(cells)
m <- rmhmodel(f)

f <- ppm(cells ~x)
m <- rmhmodel(f)

f <- ppm(cells ~1, Strauss(0.1))
m <- rmhmodel(f)

f <- ppm(cells ~1, StraussHard(r=0.1,hc=0.05))
m <- rmhmodel(f)

f <- ppm(cells ~1, Hardcore(0.07))
m <- rmhmodel(f)

f <- ppm(cells ~1, DiggleGratton(0.05,0.1))
m <- rmhmodel(f)

f <- ppm(cells ~1, Softcore(0.5), correction="isotropic")
m <- rmhmodel(f)

f <- ppm(cells ~1, Geyer(0.07,2))
m <- rmhmodel(f)

f <- ppm(cells ~1, BadGey(c(0.07,0.1,0.13),2))
m <- rmhmodel(f)

f <- ppm(cells ~1, PairPiece(r = c(0.05, 0.1, 0.2)))
m <- rmhmodel(f)

f <- ppm(cells ~1, AreaInter(r=0.06))
m <- rmhmodel(f)

# multitype

r <- matrix(0.07, 2, 2)
f <- ppm(amacrine ~1, MultiStrauss(c("off","on"),r))
m <- rmhmodel(f)

h <- matrix(min(nndist(amacrine))/2, 2, 2)
f <- ppm(amacrine ~1, MultiStraussHard(c("off","on"),r, h))
m <- rmhmodel(f)

diag(r) <- NA
diag(h) <- NA
f <- ppm(amacrine ~1, MultiStrauss(c("off","on"),r))
m <- rmhmodel(f)

f <- ppm(amacrine ~1, MultiStraussHard(c("off","on"),r, h))
m <- rmhmodel(f)

# multitype data, interaction not dependent on type

f <- ppm(amacrine ~marks, Strauss(0.05))
m <- rmhmodel(f)

# trends

f <- ppm(cells ~x, Strauss(0.1))
m <- rmhmodel(f)

f <- ppm(cells ~y, StraussHard(r=0.1,hc=0.05))
m <- rmhmodel(f)

f <- ppm(cells ~x+y, Hardcore(0.07))
m <- rmhmodel(f)

f <- ppm(cells ~polynom(x,y,2), Softcore(0.5), correction="isotropic")
m <- rmhmodel(f)

# covariates

Z <- as.im(function(x,y){ x^2+y^2 }, as.owin(cells))
f <- ppm(cells ~z, covariates=list(z=Z))
m <- rmhmodel(f)
m <- rmhmodel(f, control=list(p=1))

Zim <- as.im(Z, as.owin(cells))
f <- ppm(cells ~z, covariates=list(z=Zim))
m <- rmhmodel(f)

Z <- as.im(function(x,y){ x^2+y }, as.owin(amacrine))
f <- ppm(amacrine ~z + marks, covariates=list(z=Z))
m <- rmhmodel(f)
m <- rmhmodel(f, control=list(p=1))
m <- rmhmodel(f, control=list(p=1,fixall=TRUE))

Zim <- as.im(Z, as.owin(amacrine))
f <- ppm(amacrine ~z + marks, covariates=list(z=Zim))
m <- rmhmodel(f)

})
#
#    tests/rmhmodelHybrids.R
#
#  Test that rmhmodel.ppm and rmhmodel.default
#  work on Hybrid interaction models
#
#   $Revision: 1.4 $  $Date: 2015/12/29 08:54:49 $
#

require(spatstat)

local({
  # ......... rmhmodel.ppm .......................

  fit1 <- ppm(redwood ~1,
              Hybrid(A=Strauss(0.02), B=Geyer(0.1, 2), C=Geyer(0.15, 1)))
  m1 <- rmhmodel(fit1)
  m1
  reach(m1)

  # Test of handling 'IsOffset' 
  fit2 <- ppm(cells ~1, Hybrid(H=Hardcore(0.05), G=Geyer(0.15, 2)))
  rmhmodel(fit2)

  # Test of handling Poisson components
  fit3 <- ppm(cells ~1, Hybrid(P=Poisson(), S=Strauss(0.05)))
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
#  $Revision: 1.2 $ $Date: 2015/12/29 08:54:49 $
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
   fit <- ppm(X ~1, Poisson())
   Xsim <- rmh(fit)
   # Strauss process   
   fit <- ppm(X ~1, Strauss(r=7))
   Xsim <- rmh(fit)

   # Strauss process simulated on a larger window
   # then clipped to original window
   Xsim <- rmh(fit, control=list(nrep=Nrep, expand=1.1, periodic=TRUE))

   # Extension of model to another window (thanks to Tuomas Rajala)
   Xsim <- rmh(fit, w=square(2))
   Xsim <- simulate(fit, w=square(2))
   
   # Strauss - hard core process
#   fit <- ppm(X ~1, StraussHard(r=7,hc=2))
#   Xsim <- rmh(fit, start=list(n.start=X$n))

   # Geyer saturation process
#   fit <- ppm(X ~1, Geyer(r=7,sat=2))
#   Xsim <- rmh(fit, start=list(n.start=X$n))

   # Area-interaction process
     fit <- ppm(X ~1, AreaInter(r=7))
     Xsim <- rmh(fit, start=list(n.start=X$n))
  
     # soft core interaction process
#     X <- quadscheme(X, nd=50)
#     fit <- ppm(X ~1, Softcore(kappa=0.1), correction="isotropic")
#     Xsim <- rmh(fit, start=list(n.start=X$n))

     # Diggle-Gratton pairwise interaction model
#     fit <- ppm(cells ~1, DiggleGratton(0.05, 0.1))
#     Xsim <- rmh(fit, start=list(n.start=cells$n))
#     plot(Xsim, main="simulation from fitted Diggle-Gratton model")
   
   X <- rSSI(0.05, 100)

   # piecewise-constant pairwise interaction function
   fit <- ppm(X ~1, PairPiece(seq(0.02, 0.1, by=0.01)))
   Xsim <- rmh(fit)

   # marked point pattern
   Y <- amacrine

   # marked Poisson models
   fit <- ppm(Y)
   Ysim <- rmh(fit)

   fit <- ppm(Y~marks)
   Ysim <- rmh(fit)

   fit <- ppm(Y~x)
   Ysim <- rmh(fit)
#   fit <- ppm(Y~polynom(x,2))
#   Ysim <- rmh(fit)

   fit <- ppm(Y~marks+x)
   Ysim <- rmh(fit)
#   fit <- ppm(Y~marks+polynom(x,2))
#   Ysim <- rmh(fit)

   # multitype Strauss models
   MS <- MultiStrauss(types = levels(Y$marks),
                      radii=matrix(0.07, ncol=2, nrow=2))

#   fit <- ppm(Y~marks*polynom(x,2), MS)
    fit <- ppm(Y~marks*x, MS)
   Ysim <- rmh(fit)

   spatstat.options(op)
 })
#
#  tests/segments.R
#
#  $Revision: 1.8 $  $Date: 2015/12/29 08:54:49 $

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
d <- distppll(p, l, mintype=2, method="C")
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
#
## tests/sigtraceprogress.R
#
## Tests of *.sigtrace and *.progress
#
## $Revision: 1.3 $ $Date: 2015/12/29 08:54:49 $

require(spatstat)
local({
  plot(dclf.sigtrace(redwood, nsim=19, alternative="greater", rmin=0.02,
                     verbose=FALSE))
  plot(dclf.progress(redwood, nsim=19, alternative="greater", rmin=0.02,
                     verbose=FALSE))
  plot(dg.sigtrace(redwood, nsim=5, alternative="greater", rmin=0.02,
                     verbose=FALSE))
  plot(dg.progress(redwood, nsim=5, alternative="greater", rmin=0.02,
                   verbose=FALSE))
  ## test 'leave-two-out' algorithm
  a <- dclf.sigtrace(redwood, Lest, nsim=9, use.theory=FALSE, leaveout=2,
                     verbose=FALSE)
  aa <- dclf.progress(redwood, Lest, nsim=9, use.theory=FALSE, leaveout=2,
                      verbose=FALSE)
  b <- dg.sigtrace(redwood, Lest, nsim=5, use.theory=FALSE, leaveout=2)
  bb <- dg.progress(redwood, Lest, nsim=5, use.theory=FALSE, leaveout=2,
                    verbose=FALSE)
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


#'    tests/sparse3Darrays.R
#'  Basic tests of sparse3array.R code
#'  $Revision: 1.5 $ $Date: 2016/03/06 02:24:57 $

require(spatstat)
local({

  if(require(Matrix)) {
    M <- sparse3Darray(i=1:4, j=sample(1:4, replace=TRUE),
                       k=c(1,2,1,2), x=1:4, dims=c(5,5,2))

    M

    dimnames(M) <- list(letters[1:5], LETTERS[1:5], c("yes", "no"))
    M
    
    U <- aperm(M, c(1,3,2))
    U
    
    M[ 3:4, , ]
    
    M[ 3:4, 2:4, ]
    
    M[, 3, ]

    M[, 3, , drop=FALSE]
    
    MA <- as.array(M)
    UA <- as.array(U)

    ## tests of "[<-.sparse3Darray"
    Mflip <- Mzero <- MandM <- M
    Mflip[ , , 2:1] <- M
    stopifnot(Mflip[3,1,1] == M[3,1,2])
    Mzero[1:3,1:3,] <- 0
    stopifnot(all(Mzero[1,1,] == 0))
    M2a <- M[,,2,drop=FALSE]
    M2d <- M[,,2,drop=TRUE]
    MandM[,,1] <- M2a
    MandM[,,1] <- M2d

    ## tests of arithmetic (Math, Ops, Summary)
    negM <- -M
    oneM <- 1 * M
    twoM <- M + M
    range(M)

    cosM <- cos(M)  # non-sparse
    sinM <- sin(M)  # sparse
    
    stopifnot(all((M+M) == 2*M))     # non-sparse
    stopifnot(!any((M+M) != 2*M))    # sparse

    ## tensor operator

    tenseur(c(1,-1), M, 1, 3)
    tenseur(M, M, 1:2, 1:2)
    tenseur(M, M, 1:2, 2:1)
    V <- sparseVector(i=c(1,3,6),x=1:3, length=7)
    tenseur(V,V)
    tenseur(V,V,1,1)

    ## test of anyNA method
    anyNA(M)
    
    ## a possible application in spatstat
    cl10 <- as.data.frame(closepairs(cells, 0.1))
    cl12 <- as.data.frame(closepairs(cells, 0.12))
    cl10$k <- 1
    cl12$k <- 2
    cl <- rbind(cl10, cl12)
    n <- npoints(cells)
    Z <- with(cl,
              sparse3Darray(i=i, j=j, k=k, x=1, dims=c(n,n,2)))
    dimnames(Z) <- list(NULL, NULL, c("r=0.1", "r=0.12"))

    Z <- aperm(Z, c(3,1,2))
    stopifnot(all(sumsymouterSparse(Z) == sumsymouter(as.array(Z))))
  }
})


#
#  tests/splitpea.R
#
#  Check behaviour of split.ppp etc
#
#  Thanks to Marcelino de la Cruz
#
#  $Revision: 1.11 $  $Date: 2016/03/05 01:33:47 $
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

## More headaches with mark format
A <- runifpoint(10)
B <- runifpoint(10)
AB <- split(superimpose(A=A, B=B))

#' check that split<- respects ordering where possible
X <- amacrine
Y <- split(X)
split(X) <- Y
stopifnot(identical(X, amacrine))

#' split.ppx
df <- data.frame(x=runif(4),y=runif(4),t=runif(4),
                 age=rep(c("old", "new"), 2),
                 mineral=factor(rep(c("Au","Cu"), each=2),
                                levels=c("Au", "Cu", "Pb")),
                 size=runif(4))
X <- ppx(data=df, coord.type=c("s","s","t","m", "m","m"))
Y <- split(X, "age")
Y <- split(X, "mineral", drop=TRUE)

})
#
#   tests/step.R
#
#   $Revision: 1.4 $  $Date: 2015/12/29 08:54:49 $
#
# test for step() operation
#
require(spatstat)
local({
  Z <- as.im(function(x,y){ x^3 - y^2 }, nztrees$window)
  fitP <- ppm(nztrees ~x+y+Z, covariates=list(Z=Z))
  step(fitP)
  fitS <- update(fitP, Strauss(7))
  step(fitS)
  fitM <- ppm(amacrine ~ marks*(x+y),
              MultiStrauss(types=levels(marks(amacrine)), radii=matrix(0.04, 2, 2)))
  step(fitM)
})


##
## tests/symbolmaps.R
##
##   Quirks associated with symbolmaps, etc.
##
## $Revision: 1.3 $ $Date: 2015/12/29 08:54:49 $

local({
  require(spatstat)
  set.seed(100)
  
  ## spacing too large for tiles - upsets various pieces of code
  V <- as.im(dirichlet(runifpoint(8)))
  textureplot(V, spacing=2)

  g1 <- symbolmap(range=c(0,100), size=function(x) x/50)
  invoke.symbolmap(g1, 50, x=numeric(0), y=numeric(0), add=TRUE)

})
#
#   tests/testaddvar.R
#
# test addvar options
#
#   $Revision: 1.2 $  $Date: 2015/12/29 08:54:49 $

X <-  rpoispp(function(x,y){exp(3+3*x)})
model <- ppm(X ~y)
addvar(model, "x", crosscheck=TRUE)
addvar(model, "x", bw.input="quad")
w <- square(0.5)
addvar(model, "x", subregion=w)
addvar(model, "x", subregion=w, bw.input="points")
#
#   tests/testparres.R
#
# additional test of parres
#
#  $Revision: 1.2 $  $Date: 2015/12/29 08:54:49 $
#
require(spatstat)
local({
X <-  rpoispp(function(x,y){exp(3+x+2*x^2)})
model <- ppm(X ~x+y)

# options in parres
parres(model, "x")
parres(model, "x", bw.input="quad")
w <- square(0.5)
parres(model, "x", subregion=w)
parres(model, "x", subregion=w, bw.input="quad")

# check whether 'update.ppm' has messed up internals
mod2 <- update(model, ~x)
parres(mod2, "x")
})
#
# tests/triplets.R
#
# test code for triplet interaction
#
# $Revision: 1.5 $ $Date: 2015/12/29 08:54:49 $
#
require(spatstat)
local({
  fit <- ppm(redwood ~1, Triplets(0.1))
  fit
  suffstat(fit)
  # hard core (zero triangles, coefficient is NA)
  fit0 <- ppm(cells ~1, Triplets(0.05))
  fit0
  suffstat(fit0)
  # bug case (1 triangle in data)
  fit1 <- ppm(cells ~1, Triplets(0.15))
  fit1
  suffstat(fit1)
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




##
##  tests/updateppm.R
##
##  Check validity of update.ppm
##
##  $Revision: 1.4 $ $Date: 2016/03/08 06:30:46 $

local({
    require(spatstat)
    h <- function(m1, m2) {
        mc <- deparse(sys.call())
        cat(paste(mc, "\t... "))
        m1name <- deparse(substitute(m1))
        m2name <- deparse(substitute(m2))
        if(!identical(names(coef(m1)), names(coef(m2))))
            stop(paste("Differing results for", m1name, "and", m2name,
                       "in updateppm.R"),
                 call.=FALSE)
        cat("OK\n")
    }
    X <- redwood[c(TRUE,FALSE)]
    Y <- redwood[c(FALSE,TRUE)]
    fit0f <- ppm(X ~ 1, nd=8)
    fit0p <- ppm(X, ~1, nd=8)
    fitxf <- ppm(X ~ x, nd=8)
    fitxp <- ppm(X, ~x, nd=8)

    cat("Basic consistency ...\n")
    h(fit0f, fit0p)
    h(fitxf, fitxp)

    cat("\nTest correct handling of model formulas ...\n")
    h(update(fitxf, Y), fitxf)
    h(update(fitxf, Q=Y), fitxf)
    h(update(fitxf, Y~x), fitxf)
    h(update(fitxf, Q=Y~x), fitxf)
    h(update(fitxf, ~x), fitxf)

    h(update(fitxf, Y~1), fit0f)
    h(update(fitxf, ~1), fit0f)
    h(update(fit0f, Y~x), fitxf)
    h(update(fit0f, ~x), fitxf)

    h(update(fitxp, Y), fitxp)
    h(update(fitxp, Q=Y), fitxp)
    h(update(fitxp, Y~x), fitxp)
    h(update(fitxp, Q=Y~x), fitxp)
    h(update(fitxp, ~x), fitxp)

    h(update(fitxp, Y~1), fit0p)
    h(update(fitxp, ~1), fit0p)
    h(update(fit0p, Y~x), fitxp)
    h(update(fit0p, ~x), fitxp)

    cat("\nTest scope handling for left hand side ...\n")
    X <- Y
    h(update(fitxf), fitxf)

    cat("\nTest scope handling for right hand side ...\n")
    Z <- distmap(X)
    fitZf <- ppm(X ~ Z)
    fitZp <- ppm(X, ~ Z)
    h(update(fitxf, X ~ Z), fitZf)
    h(update(fitxp, X ~ Z), fitZp)
    h(update(fitxf, . ~ Z), fitZf)
    h(update(fitZf, . ~ x), fitxf)
    h(update(fitZf, . ~ . - Z), fit0f)
    h(update(fitxp, . ~ Z), fitZp)
    h(update(fitZp, . ~ . - Z), fit0p)
    h(update(fit0p, . ~ . + Z), fitZp)
    h(update(fitZf, . ~ . ), fitZf)
    h(update(fitZp, . ~ . ), fitZp)

    cat("\nTest use of internal data ...\n")
    h(update(fitZf, ~ x, use.internal=TRUE), fitxf)
    fitsin <- update(fitZf, X~sin(Z))
    h(update(fitZf, ~ sin(Z), use.internal=TRUE), fitsin)

    cat("\nTest step() ... ")
    fut <- ppm(X ~ Z + x + y, nd=8)
    fut0 <- step(fut, trace=0)
    cat("OK\n")
})

# test update.lppm

local({
  X <- runiflpp(20, simplenet)
  fit0 <- lppm(X ~ 1)
  fit1 <- update(fit0, ~ x)
  anova(fit0, fit1, test="LR")
  cat("update.lppm(fit, ~trend) is OK\n")
  fit2 <- update(fit0, . ~ x)
  anova(fit0, fit2, test="LR")
  cat("update.lppm(fit, . ~ trend) is OK\n")
})
#
#  tests/vcovppm.R
#
#  Check validity of vcov.ppm algorithms
#
#  Thanks to Ege Rubak
#
#  $Revision: 1.6 $  $Date: 2015/12/29 08:54:49 $
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
  modelY <- ppm(Y ~1, Geyer(0.1, 1))

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


  ## tests of 'deltasuffstat' code
  ##     Handling of offset terms
  modelH <- ppm(cells ~x, Hardcore(0.05))
  a <- vcov(modelH, generic=TRUE) ## may fall over
  b <- vcov(modelH, generic=FALSE)
  if(disagree(a, b))
    stop("Disagreement between vcov.ppm algorithms for Hardcore model")
  
  ##     Correctness of pairwise.family$delta2
  modelZ <- ppm(amacrine ~1, MultiStrauss(radii=matrix(0.1, 2, 2)))
  b <- vcov(modelZ, generic=FALSE)
  g <- vcov(modelZ, generic=TRUE)
  if(disagree(b, g))
    stop("Disagreement between vcov.ppm algorithms for MultiStrauss model")

  ## Test that 'deltasuffstat' works for Hybrids
  modHyb <- ppm(japanesepines ~ 1, Hybrid(Strauss(0.05), Strauss(0.1)))
})
#
# tests/windows.R
#
# Tests of owin geometry code
#
#  $Revision: 1.3 $  $Date: 2015/12/29 08:54:49 $

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

  ## thanks to Brian Ripley / SpatialVx
  M <- as.mask(letterR)
  stopifnot(area(bdry.mask(M)) > 0)
  stopifnot(area(convexhull(M)) > 0)
  R <- as.mask(square(1))
  stopifnot(area(bdry.mask(R)) > 0)
  stopifnot(area(convexhull(R)) > 0)
})

##
## tests/xysegment.R
##
##    Test weird problems and boundary cases for line segment code
##
##    $Version$ $Date: 2016/02/12 08:18:08 $ 
##
require(spatstat)
local({
  # segment of length zero
  B <- psp(1/2, 1/2, 1/2, 1/2, window=square(1))
  BB <- angles.psp(B)
  A <- runifpoint(3)
  AB <- project2segment(A,B)

  # mark inheritance
  X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  marks(X) <- 1:10
  Y <- selfcut.psp(X)
  marks(X) <- data.frame(A=1:10, B=factor(letters[1:10]))
  Z <- selfcut.psp(X)
})
