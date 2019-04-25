#'
#'    tests/quadschemes.R
#'
#'    Quadrature schemes, dummy points etc
#' 
#'   $Revision: 1.6 $ $Date: 2019/01/20 05:49:40 $
#'

require(spatstat)
local({
  ##  class 'quad' 
  qu <- quadscheme(cells)
  qm <- quadscheme(amacrine)
  plot(qu)
  plot(qm)
  is.multitype(qu)
  is.multitype(qm)
  a <- param.quad(qu)
  a <- param.quad(qm)
  a <- equals.quad(qu)
  a <- equals.quad(qm)
  a <- domain(qu)
  unitname(qu) <- c("Furlong", "Furlongs")
  
  ## utilities
  b <- cellmiddles(square(1), 3, 4)
  b <- cellmiddles(letterR, 3, 4, distances=FALSE)
  b <- cellmiddles(letterR, 3, 4, distances=TRUE)
  v <- tilecentroids(square(1), 3, 4)
  v <- tilecentroids(letterR, 3, 4)
  n <- default.n.tiling(cells)
  n <- default.n.tiling(cells, nd=4)
  n <- default.n.tiling(cells, ntile=4)
  n <- default.n.tiling(cells, ntile=4, quasi=TRUE)

  ## quadrature weights - special cases
  X <- runifpoint(10, as.mask(letterR))
  gr <- gridweights(X, ntile=12, npix=7) # causes warnings about zero digital area
  
  ## plot.quad 
  plot(quadscheme(cells, method="dirichlet", nd=7),              tiles=TRUE)
  plot(quadscheme(cells, method="dirichlet", nd=7, exact=FALSE), tiles=TRUE)
  
  ## logistic
  d <- quadscheme.logi(cells, logi.dummy(cells, "binomial"))
  print(summary(d))
  d <- quadscheme.logi(cells, logi.dummy(cells, "poisson"))
  print(summary(d))
  d <- quadscheme.logi(cells, logi.dummy(cells, "grid"))
  print(summary(d))
  d <- quadscheme.logi(cells, logi.dummy(cells, "transgrid"))
  print(summary(d))
  d <- quadscheme.logi(amacrine,
                       logi.dummy(amacrine, "binomial", mark.repeat=TRUE))
  print(summary(d))
  d <- quadscheme.logi(amacrine,
                       logi.dummy(amacrine, "poisson", mark.repeat=FALSE))
  print(summary(d))
})
#'  tests/randoms.R
#'   Further tests of random generation code
#'  $Revision: 1.7 $ $Date: 2019/04/13 14:39:48 $

require(spatstat)
local({
  A <- runifrect(6, nsim=2)
  A <- runifdisc(6, nsim=2)
  A <- runifpoispp(5, nsim=2)
  A <- runifpoispp(0, nsim=2)
  Z <- as.im(function(x,y) 10*x, square(1))
  A <- rpoint(n=6, f=Z, fmax=10, nsim=2)
  A <- rSSI(0.05, 6, nsim=2)
  A <- rSSI(0.05, 10, win=square(c(-0.5, 1.5)), x.init=A[[1]], nsim=2)  
  A <- rstrat(nx=4, nsim=2)
  A <- rsyst(nx=4, nsim=2)
  A <- rthin(cells, P=0.5, nsim=2)
  A <- rthin(cells, runif(42))
  A <- rjitter(cells, nsim=2, retry=FALSE)

  op <- spatstat.options(fastpois=FALSE)
  A <- runifpoispp(5, nsim=2)
  A <- rpoispp(Z)
  spatstat.options(op)

  b3 <- box3(c(0,1))
  b4 <- boxx(c(0,1), c(0,1), c(0,1), c(0,1))
  X <- rMaternInhibition(2, kappa=20, r=0.1, win=b3)
  Y <- rMaternInhibition(2, kappa=20, r=0.1, win=b4)
})

local({
  f1 <- function(x,y){(x^2 + y^3)/10}
  f2 <- function(x,y){(x^3 + y^2)/10}
  ZZ <- solist(A=as.im(f1, letterR),
               B=as.im(f2, letterR))
  XX <- rmpoispp(ZZ, nsim=3)
  YY <- rmpoint(10, f=ZZ, nsim=3)
  g <- function(x,y,m){(10+as.integer(m)) * (x^2 + y^3)}
  VV <- rpoint.multi(10, f=g,
                     marks=factor(sample(letters[1:3], 10, replace=TRUE)),
                     nsim=3)
  L <- edges(letterR)
  E <- runifpoisppOnLines(5, L)
  G <- rpoisppOnLines(ZZ, L)
})

local({
  #' cluster models + bells + whistles
  X <- rThomas(10, 0.2, 5, saveLambda=TRUE)
  X <- rMatClust(10, 0.05, 4, saveLambda=TRUE)
  X <- rCauchy(30, 0.01, 5, saveLambda=TRUE)
  X <- rVarGamma(30, 2, 0.02, 5, saveLambda=TRUE)
  Z <- as.im(function(x,y){ 5 * exp(2 * x - 1) }, owin())
  Y <- rThomas(10, 0.2, Z, saveLambda=TRUE)
  Y <- rMatClust(10, 0.05, Z, saveLambda=TRUE)
  Y <- rCauchy(30, 0.01, Z, saveLambda=TRUE)
  Y <- rVarGamma(30, 2, 0.02, Z, saveLambda=TRUE)
})

reset.spatstat.options()


#'  tests/resid.R
#'
#'  Stuff related to residuals and residual diagnostics
#'
#'   $Revision: 1.2 $  $Date: 2018/06/11 07:07:31 $
#'

require(spatstat)
local({
  fit <- ppm(cells ~x, Strauss(r=0.15))
  diagnose.ppm(fit, cumulative=FALSE)
  diagnose.ppm(fit, cumulative=FALSE, type="pearson")

  fitoff <- ppm(cells ~ sin(x) + offset(y))
  plot(a <- parres(fitoff, "x"))
  plot(b <- parres(fitoff, "y"))
  print(a)
  print(b)

  d <- diagnose.ppm(fit, which="marks")
  plot(d, plot.neg="discrete")
  plot(d, plot.neg="imagecontour")

  d <- diagnose.ppm(fit, type="pearson", which="smooth")
  plot(d, plot.smooth="image")
  plot(d, plot.smooth="contour")
  plot(d, plot.smooth="imagecontour")
  
  d <- diagnose.ppm(fit, type="pearson", which="x")
  plot(d)
  d <- diagnose.ppm(fit, type="pearson", which="y")
  plot(d)
  
  diagnose.ppm(fit, type="pearson", which="x", cumulative=FALSE)
  diagnose.ppm(fit, type="pearson", which="x", cumulative=FALSE)
  diagnose.ppm(fit, type="raw", plot.neg="discrete", plot.smooth="image")
  diagnose.ppm(fit, type="pearson", plot.neg="contour", plot.smooth="contour")

  diagnose.ppm(fitoff, type="raw", which="smooth", plot.smooth="persp")
  diagnose.ppm(fitoff, type="pearson", plot.neg="imagecontour")

  plot(Frame(letterR), main="")
  ploterodewin(letterR, erosion(letterR, 0.05), main="jeans")
  W <- as.mask(letterR)
  plot(Frame(W), main="")
  ploterodewin(W, erosion(W, 0.05), main="JeAnS")
})


##
## tests/rhohat.R
##
## Test all combinations of options for rhohatCalc
##
## $Revision: 1.3 $ $Date: 2018/05/13 04:42:21 $

local({
  require(spatstat)
  X <-  rpoispp(function(x,y){exp(3+3*x)})
  ## rhohat.ppp
  ## done in example(rhohat):
  ## rhoA <- rhohat(X, "x")
  ## rhoB <- rhohat(X, "x", method="reweight")
  ## rhoC <- rhohat(X, "x", method="transform")

  ## alternative smoother (if package locfit available)
  rhoA <- rhohat(X, "x", smoother="local")
  rhoB <- rhohat(X, "x", smoother="local", method="reweight")
  rhoC <- rhohat(X, "x", smoother="local", method="transform")

  ## rhohat.ppm
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

  ## class support
  plot(rhoA)
  plot(rhoA, rho ~ x, shade=NULL)
  plot(rhoA, log(rho) ~ x, shade=NULL)
  plot(rhoA, log(.) ~ x)

  ## rho2hat
  r2xy <- rho2hat(X, "x", "y")
  r2xyw <- rho2hat(X, "x", "y", method="reweight")
  plot(r2xy, do.points=TRUE)
  xcoord <- function(x,y) x
  ycoord <- function(x,y) y
  xim <- as.im(xcoord, W=Window(X))
  r2fi <- rho2hat(X, ycoord, xim)
  r2if <- rho2hat(X, xim, ycoord)
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
                   nburn=0, nsave=1, pstage="start"))
   #' Need to set pstage='start' so that proposals are generated
   #' at the start of the procedure in both cases.

   stopifnot(npoints(X.save) == npoints(X.nosave))
   stopifnot(max(nncross(X.save, X.nosave)$dist) == 0)
   stopifnot(max(nncross(X.nosave, X.save)$dist) == 0)
})
##
##   tests/rmhBasic.R
##
##   $Revision: 1.17 $  $Date: 2019/04/25 03:11:55 $
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
   X1.strauss2 <- rmh(model=mod01,start=list(n.start=80),
                     control=list(nrep=nr, periodic=FALSE))

   # Strauss process, conditioning on n = 80:
   X2.strauss <- rmh(model=mod01,start=list(n.start=80),
                     control=list(p=1,nrep=nr))
   stopifnot(npoints(X2.strauss) == 80)

   # test tracking mechanism
   X1.strauss <- rmh(model=mod01,start=list(n.start=80),
                     control=list(nrep=nr), track=TRUE)
   X2.strauss <- rmh(model=mod01,start=list(n.start=80),
                         control=list(p=1,nrep=nr), track=TRUE)
   
   # Hard core process:
   mod02 <- list(cif="hardcore",par=list(beta=2,hc=0.7),w=c(0,10,0,10))
   X3.hardcore <- rmh(model=mod02,start=list(n.start=60),
                     control=list(nrep=nr))
   X3.hardcore2 <- rmh(model=mod02,start=list(n.start=60),
                       control=list(nrep=nr, periodic=FALSE))

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
   stopifnot(npoints(X5.strauss) == 42)

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
   X1.straush2 <- rmh(model=mod04,start=list(n.start=70),
                     control=list(nrep=nr, periodic=FALSE))
   
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

   # Fiksel
   modFik <- list(cif="fiksel",
                  par=list(beta=180,r=0.15,hc=0.07,kappa=2,a= -1.0),
                  w=square(1))
   X.fiksel <- rmh(model=modFik,start=list(n.start=10),
                   control=list(nrep=nr))
   X.fiksel2 <- rmh(model=modFik,start=list(n.start=10),
                   control=list(nrep=nr,periodic=FALSE))

   # Penttinen process:
   modpen <- rmhmodel(cif="penttinen",par=list(beta=2,gamma=0.6,r=1),
                 w=c(0,10,0,10))
   X.pen <- rmh(model=modpen,start=list(n.start=10),
                control=list(nrep=nr))
   X.pen2 <- rmh(model=modpen,start=list(n.start=10),
                 control=list(nrep=nr, periodic=FALSE))

   # Area-interaction, inhibitory
   mod.area <- list(cif="areaint",par=list(beta=2,eta=0.5,r=0.5), w=square(10))
   X.area <- rmh(model=mod.area,start=list(n.start=60),
                 control=list(nrep=nr))
   X.areaE <- rmh(model=mod.area,start=list(n.start=60),
                 control=list(nrep=nr, periodic=FALSE))

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
   X.sftcr2 <- rmh(model=mod07,start=list(n.start=70),
                   control=list(nrep=nr, periodic=FALSE))
   
   # Diggle, Gates, and Stibbard:
   mod12 <- list(cif="dgs",par=list(beta=3600,rho=0.08),w=c(0,1,0,1))
   X.dgs <- rmh(model=mod12,start=list(n.start=300),
                control=list(nrep=nr))
   X.dgs2 <- rmh(model=mod12,start=list(n.start=300),
                 control=list(nrep=nr, periodic=FALSE))
   
   # Diggle-Gratton:
   mod13 <- list(cif="diggra",
                 par=list(beta=1800,kappa=3,delta=0.02,rho=0.04),
                 w=square(1))
   X.diggra <- rmh(model=mod13,start=list(n.start=300),
                   control=list(nrep=nr))
   X.diggra2 <- rmh(model=mod13,start=list(n.start=300),
                    control=list(nrep=nr, periodic=FALSE))
   
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
   X2.geyer2 <- rmh(model=mod15,start=list(n.start=200),
                   control=list(nrep=nr, periodic=FALSE))
   
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
                      control=list(nrep=nr, periodic=TRUE))
      X.lookup2 <- rmh(model=mod17,start=list(n.start=100),
                       control=list(nrep=nr, expand=1, periodic=FALSE))

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

   #' Test other code blocks
   #'  nsim > 1
   Xlist <- rmh(model=mod01,start=list(n.start=80),
             control=list(nrep=nr),
             nsim=2)
   #' Condition on contents of window
   XX <- Xlist[[1]]
   YY <- XX[square(2)]
   XXwindow <- rmh(model=mod01, start=list(n.start=80),
                   control=list(nrep=nr, x.cond=YY))
   #' Palm conditioning
   XXpalm <- rmh(model=mod01,start=list(n.start=80),
             control=list(nrep=nr, x.cond=coords(YY)))


   #' code blocks for various interactions, not otherwise tested
   rr <- seq(0,0.2,length=8)[-1]
   gmma <- c(0.5,0.6,0.7,0.8,0.7,0.6,0.5)
   mod18 <- list(cif="badgey",par=list(beta=4000, gamma=gmma,r=rr,sat=5),
                 w=square(1))
   Xbg <- rmh(model=mod18,start=list(n.start=20),
              control=list(nrep=1e4, periodic=TRUE))
   Xbg2 <- rmh(model=mod18,start=list(n.start=20),
              control=list(nrep=1e4, periodic=FALSE))

})

reset.spatstat.options()
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
#  $Revision: 1.4 $  $Date: 2018/07/21 00:46:16 $
#

require(spatstat)
local({
fit <- ppm(cells ~x)

# check rmhmodel.ppm
mod <- rmhmodel(fit)
is.expandable(mod)
wsim <- as.rectangle(mod$trend)
# work around changes in 'unitname'
wcel <- as.owin(cells)
unitname(wcel) <- unitname(cells)
# test
if(!identical(wsim, wcel))
  stop("Expansion occurred improperly in rmhmodel.ppm")
})


#
#  tests/rmhMulti.R
#
#  tests of rmh, running multitype point processes
#
#   $Revision: 1.8 $  $Date: 2019/04/17 09:12:41 $

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

   # Multitype hardcore:
   rhc  <- matrix(c(9.1,5.0,5.0,2.5),2,2)
   mod085 <- list(cif="multihard",par=list(beta=beta,hradii=rhc),
                  w=c(0,250,0,250))
   X.multihard <- rmh(model=mod085,start=list(n.start=80),
                      control=list(ptypes=c(0.75,0.25),nrep=nr,nverb=nv))
   X.multihardP <- rmh(model=mod085,start=list(n.start=80),
                       control=list(ptypes=c(0.75,0.25),nrep=nr,nverb=nv,
                                    periodic=TRUE))
   #' check handling of initial state which violates hard core
   Xclose <- ppp(c(10,12,20),c(10,10,10), marks=factor(c(1,1,2)),
              c(0,250),c(0,250))
   X.multihard.close <- rmh(model=mod085,start=list(x.start=Xclose),
                            control=list(ptypes=c(0.75,0.25),nrep=nr,nverb=nv))
   X.multihard.closeP <- rmh(model=mod085,start=list(x.start=Xclose),
                             control=list(ptypes=c(0.75,0.25),nrep=nr,nverb=nv,
                                          periodic=TRUE))

   # Multitype Strauss hardcore:
   mod09 <- list(cif="straushm",par=list(beta=beta,gamma=gmma,
                iradii=r,hradii=rhc),w=c(0,250,0,250))
   X.straushm <- rmh(model=mod09,start=list(n.start=80),
                     control=list(ptypes=c(0.75,0.25),nrep=nr,nverb=nv))
   X.straushm.close <- rmh(model=mod09,start=list(x.start=Xclose),
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

reset.spatstat.options()
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

  ## Test of handling 'IsOffset' 
  fit2 <- ppm(cells ~1, Hybrid(H=Hardcore(0.05), G=Geyer(0.15, 2)))
  m2 <- rmhmodel(fit2)
  ## also test C code for hybrid interaction with hard core
  fakecells <- rmh(fit2, nrep=1e4)

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
#  $Revision: 1.4 $ $Date: 2019/02/21 01:59:48 $
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
  
   # Penttinen process
     fit <- ppm(X ~1, Penttinen(r=7))
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

   #' marked Poisson models
   fit <- ppm(Y)
   Ysim <- rmh(fit)

   fit <- ppm(Y~marks)
   Ysim <- rmh(fit)

   fit <- ppm(Y~x)
   Ysim <- rmh(fit)

   fit <- ppm(Y~marks+x)
   Ysim <- rmh(fit)

   #' multitype Strauss
   typ <- levels(Y$marks)
   MS <- MultiStrauss(types = typ,
                      radii=matrix(0.07, ncol=2, nrow=2))

   fit <- ppm(Y~marks*x, MS)
   Ysim <- rmh(fit)

   #' multitype Hardcore
   h0 <- minnndist(unmark(Y)) * 0.95
   MH <- MultiHard(types = typ,
                   hradii=matrix(h0, ncol=2, nrow=2))
   fit <- ppm(Y ~ marks+x, MH)
   Ysim <- rmh(fit)
   #' other code blocks
   Ysim <- rmh(fit, control=list(periodic=TRUE, expand=1))
   Ysim <- rmh(fit, control=list(periodic=FALSE, expand=1))
   #' multihard core with invalid initial state
   Ydouble <- superimpose(Y, rjitter(Y, h0/10))
   Ysim <- rmh(fit, start=list(x.start=Ydouble))

   #' Lennard-Jones
   fut <- ppm(unmark(longleaf) ~ 1, LennardJones(), rbord=1)
   Ysim <- rmh(fut)
   Ysim <- rmh(fut, control=list(periodic=TRUE, expand=1))
   
   spatstat.options(op)
 })


reset.spatstat.options()
#'
#'   tests/rmhsnoopy.R
#'
#'   Test the rmh interactive debugger
#' 
#'   $Revision: 1.8 $  $Date: 2018/10/17 08:57:32 $

require(spatstat)
local({
  ## fit a model and prepare to simulate
  R <- 0.1
  fit <- ppm(cells ~ 1, Strauss(R))
  siminfo <- rmh(fit, preponly=TRUE)
  Wsim <- siminfo$control$internal$w.sim
  Wclip <- siminfo$control$internal$w.clip
  if(is.null(Wclip)) Wclip <- Window(cells)

  ## determine debugger interface panel geometry
  P <- rmhsnoop(Wsim=Wsim, Wclip=Wclip, R=R,
                xcoords=runif(40),
                ycoords=runif(40),
                mlevels=NULL, mcodes=NULL,
                irep=3, itype=1,
                proptype=1, proplocn=c(0.5, 0.5), propmark=0, propindx=0,
                numerator=42, denominator=24,
                panel.only=TRUE)
  boxes <- P$boxes
  clicknames <- names(P$clicks)
  boxcentres <- do.call(concatxy, lapply(boxes, centroid.owin))

  ## design a sequence of clicks
  actionsequence <- c("Up", "Down", "Left", "Right",
                      "At Proposal", "Zoom Out", "Zoom In", "Reset",
                      "Accept", "Reject", "Print Info",
                      "Next Iteration", "Next Shift", "Next Death",
                      "Skip 10", "Skip 100", "Skip 1000", "Skip 10,000",
                      "Skip 100,000", "Exit Debugger")
  actionsequence <- match(actionsequence, clicknames)
  actionsequence <- actionsequence[!is.na(actionsequence)]
  xy <- lapply(boxcentres, "[", actionsequence)

  ## queue the click sequence
  spatstat.utils::queueSpatstatLocator(xy$x,xy$y)

  ## go
  rmh(fit, snoop=TRUE)
})
