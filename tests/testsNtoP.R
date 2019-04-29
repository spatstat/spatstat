#
# tests/NAinCov.R
#
# Testing the response to the presence of NA's in covariates
#
# $Revision: 1.5 $ $Date: 2015/12/29 08:54:49 $

require(spatstat)
local({
  X <- runifpoint(42)
  Y <- as.im(function(x,y) { x+y }, owin())
  Y[owin(c(0.2,0.4),c(0.2,0.4))] <- NA
  # fit model: should produce a warning but no failure
  misfit <- ppm(X ~Y, covariates=list(Y=Y))
  # prediction 
  Z <- predict(misfit, type="trend", se=TRUE)
  # covariance matrix: all should be silent
  v <- vcov(misfit)
  ss <- vcov(misfit, what="internals")
  NULL
})





#
#    tests/nndist.R
#
# Check that nndist and nnwhich give
# results consistent with direct calculation from pairdist
#
# Similarly for nncross and distfun
#
# Also test whether minnndist(X) == min(nndist(X))
#
# Also test nnorient()
#
#   $Revision: 1.24 $  $Date: 2019/04/29 06:34:29 $
#

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

  Y <- runifpoint3(17)
  a <- nncross(X,Y)
  a <- nncross(X,Y, what="dist")
  a <- nncross(X,Y, what="which")
  a2 <- nncross(X,Y, k=2)
  a2 <- nncross(X,Y, what="dist", k=2)
  a2 <- nncross(X,Y, what="which", k=2)
  iX <- 1:42
  iZ <- 30:42
  Z <- X[iZ]
  b <- nncross(X, Z, iX=iX, iY=iZ)
  b <- nncross(X, Z, iX=iX, iY=iZ, what="which")
  b <- nncross(X, Z, iX=iX, iY=iZ, what="dist")
  b2 <- nncross(X, Z, iX=iX, iY=iZ, k=2)
  b2 <- nncross(X, Z, iX=iX, iY=iZ, what="which", k=2)
  b2 <- nncross(X, Z, iX=iX, iY=iZ, what="dist", k=2)
  
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
  ndw <- nncross(X, Y, k=1:4, what="which")
  if(any(is.na(ndw)))
    stop("NA's returned by nncross.ppp(k > 1, what='which')")
  nnc4 <- nncross(X, Y, k=1:4)
  iswhich <- (substr(colnames(nnc4), 1, nchar("which")) == "which")
  ndw <- nnc4[,iswhich]
  if(any(is.na(ndw)))
    stop("NA's returned by nncross.ppp(k > 1)$which")
  
  # test of correctness for nncross with k > 1
  flipcells <- flipxy(cells)
  calcwhich <- nncross(cells, flipcells, k=1:4, what="which")
  truewhich <- t(apply(crossdist(cells,flipcells), 1, order))[,1:4]
  if(any(calcwhich != truewhich))
    stop("nncross(k > 1) gives wrong answer")

  #' cover some C code blocks
  Z <- runifpoint(50)
  X <- Z[1:30]
  Y <- Z[20:50]
  iX <- 1:30
  iY <- 20:50
  Ndw <- nncross(X,Y, iX, iY, k=3)
  Nw  <- nncross(X,Y, iX, iY, k=3, what="which")
  Nd  <- nncross(X,Y, iX, iY, k=3, what="dist")
  
  # test of agreement between nngrid.h and knngrid.h
  #    dimyx=23 (found by trial-and-error) ensures that there are no ties 
  a <- as.matrix(nnmap(cells, what="which", dimyx=23))
  b <- as.matrix(nnmap(cells, what="which", dimyx=23, k=1:2)[[1]])
  if(any(a != b))
    stop("algorithms in nngrid.h and knngrid.h disagree")

  ## minnndist
  mfast <- minnndist(X)
  mslow <- min(nndist(X))
  if(abs(mfast-mslow) > eps)
    stop("minnndist(X) disagrees with min(nndist(X))")
  mfast <- maxnndist(X)
  mslow <- max(nndist(X))
  if(abs(mfast-mslow) > eps)
    stop("maxnndist(X) disagrees with max(nndist(X))")
})

local({
  # tests for has.close()
  # (the default method uses nndist or pairdist, and can be trusted!)
  a <- has.close(redwood, 0.05)
  b <- has.close.default(redwood, 0.05)
  if(any(a != b)) stop("Incorrect result for has.close(X, r)")

  a <- has.close(redwood, 0.05, periodic=TRUE)
  a <- has.close.default(redwood, 0.05, periodic=TRUE)
  if(any(a != b)) stop("Incorrect result for has.close(X, r, periodic=TRUE)")

  Y <- split(amacrine)
  a <- with(Y, has.close(on, 0.05, off))
  b <- with(Y, has.close.default(on, 0.05, off))
  if(any(a != b)) stop("Incorrect result for has.close(X, r, Y)")

  a <- with(Y, has.close(on, 0.05, off, periodic=TRUE))
  b <- with(Y, has.close.default(on, 0.05, off, periodic=TRUE))
  if(any(a != b)) stop("Incorrect result for has.close(X, r, Y, periodic=TRUE)")
})

local({
  b <- bdist.pixels(letterR, style="coords")
})

local({
  #' test nnorient
  nnorient(cells, domain=erosion(Window(cells), 0.1))
  #' degenerate case
  X <- cells[nndist(cells) > bdist.points(cells)]
  f <- nnorient(X)
  #' nnclean
  A <- nnclean(shapley, k=17, edge.correct=TRUE)
  B <- nnclean(runifpoint3(300), 3)
  #' stienen set
  #' bug when disc radius is zero
  Y <- unmark(humberside)[40:100] # contains duplicated points
  stienen(Y)
  Z <- stienenSet(Y)
})

  
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

#'   tests/perspim.R
#'
#'   Check persp.im handling of NA, etc
#' 
#'   $Revision: 1.1 $  $Date: 2016/08/27 02:53:35 $

require(spatstat)

local({
  set.seed(42)
  Z <- distmap(letterR, invert=TRUE)[letterR, drop=FALSE]
  X <- runifpoint(100, Frame(Z))
  M <- persp(Z, colin=Z, visible=TRUE)
  perspPoints(X, Z=Z, M=M)
})
##
## tests/pixelgripes.R
##     Problems related to pixellation of windows
##
## $Revision: 1.4 $ $Date: 2018/10/10 08:04:10 $

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

local({
  ## pixellate.ppp includes mapping from (x,y) to (row, col)
  Z <- pixellate(cells, savemap=TRUE)
  ind <- attr(Z, "map")
  m <- (as.matrix(Z))[ind]
  if(!all(m == 1)) stop("Coordinate mismatch in pixellate.ppp")
})

## 
## tests/polygons.R
##
##  $Revision: 1.3 $ $Date: 2018/07/22 02:10:07 $
##
require(spatstat)
local({
  co <- as.ppp(corners(letterR), letterR, check=FALSE)
  co[letterR]

  b <- letterR$bdry
  a <- sapply(b, xypolyselfint, yesorno=TRUE)
  a <- lapply(b, xypolyselfint, proper=TRUE)
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

#'
#'  tests/ppmlogi.R
#'
#' Tests of ppm(method='logi')
#'    and related code (predict, leverage etc)
#'
#' $Revision: 1.11 $  $Date: 2019/04/12 03:35:12 $
#'

require(spatstat)
local({
  fit <- ppm(cells ~x, method="logi")
  f <- fitted(fit)
  p <- predict(fit)
  u <- summary(fit)
  fitS <- ppm(cells ~x, Strauss(0.12), method="logi")
  fS <- fitted(fitS)
  pS <- predict(fitS)
  uS <- summary(fitS)

  plot(leverage(fit))
  plot(influence(fit))
  plot(dfbetas(fit))
  plot(leverage(fitS))
  plot(influence(fitS))
  plot(dfbetas(fitS))
})

local({
  #' same with hard core - A1 is singular
  fitH <- ppm(cells ~x, Strauss(0.08), method="logi")
  fH <- fitted(fitH)
  pH <- predict(fitH)
  uH <- summary(fitH)
  plot(leverage(fitH))
  plot(influence(fitH))
  plot(dfbetas(fitH))
})
  
local({
  #' logistic fit to data frame of covariates
  z <- c(rep(TRUE, 5), rep(FALSE, 5))
  df <- data.frame(A=z + 2* runif(10),
                   B=runif(10))
  Y <- quadscheme.logi(runifpoint(5), runifpoint(5))
  fut <- ppm(Y ~ A+B, data=df, method="logi")
  sf <- summary(fut)
})

local({
  #' vblogit code, just to check that it runs.
  fee <- ppm(cells ~ x, method="VBlogi", nd=21)
  print(fee)
  summary(fee)
  Z <- predict(fee)
  summary(Z)
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

  ## (4) anova.ppm
  fut1  <- ppm(cells ~ 1, Strauss(0.1))
  futx  <- ppm(cells ~ x, Strauss(0.1))
  anova(fut1, test="Chi")
  anova(futx, test="Chi")
  fut1a <- ppm(cells ~ 1, Strauss(0.1), rbord=0)
  anova(fut1a, futx, test="Chi")
  fut1d <- ppm(cells ~ 1, Strauss(0.1), nd=23)
  anova(fut1d, futx, test="Chi")
  ## The following doesn't work yet
  ## futxyg <- ppm(cells ~ x + s(y), Strauss(0.1), use.gam=TRUE)
  ## anova(futx, futxyg)
  fatP <- ppm(amacrine ~ marks)
  fatM <- ppm(amacrine ~ marks, MultiStrauss(matrix(0.07, 2, 2)))
  anova(fatP, fatM, test="Chi")
})

grep#
#   tests/ppmtricks.R
#
#   Test backdoor exits and hidden options in ppm
#        and summary.ppm, print.summary.ppm
#
#   Plus assorted tricks
#
#   $Revision: 1.10 $  $Date: 2019/02/10 07:00:38 $
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
  mop <- spatstat.options(maxmatrix=5000)
  pr <- predict(ppm(cells ~ x, AreaInter(0.05)))
  spatstat.options(mop)
  
  ## (5) shortcuts in summary.ppm
  ## and corresponding behaviour of print.summary.ppm
  print(summary(fit, quick=TRUE))
  print(summary(fit, quick="entries"))
  print(summary(fit, quick="no prediction"))
  print(summary(fit, quick="no variances"))

  ## (6) suffstat.R
  fitP <- update(fit, Poisson())
  suffstat.poisson(fitP, cells)
  fit0 <- killinteraction(fit)
  suffstat.poisson(fit0, cells)

  ## (7) support for class ppm
  Z <- as.im(function(x,y){x}, Window(cells))
  fitZ <- ppm(cells ~ Z)
  U <- getppmOriginalCovariates(fitZ)
  logLik(fitZ, absolute=TRUE)
  fitNot <- ppm(redwood ~1, Strauss(0.1))
  fitFast <- emend(fitNot, trace=TRUE)
  op <- spatstat.options(project.fast=TRUE)
  fitFast <- emend(fitNot, trace=TRUE)
  spatstat.options(op)
  
  fut <- kppm(redwood ~ x)
  A <- quad.ppm(fut)

  ## (8) support for class profilepl
  rr <- data.frame(r=seq(0.05, 0.15, by=0.02))
  ps <- profilepl(rr, Strauss, cells)
  plot(ps)
  simulate(ps, nrep=1e4)
  parameters(ps)
  fitin(ps)
  predict(ps, type="cif")

  ## (9) class 'plotppm'
  fut <- ppm(amacrine ~ marks + polynom(x,y,2), Strauss(0.07))
  p <- plot(fut, plot.it=FALSE)
  print(p)
  plot(p, how="contour")
  plot(p, how="persp")
})

reset.spatstat.options()
#'
#'   tests/ppp.R
#'
#'   $Revision: 1.4 $ $Date: 2018/11/05 00:50:01 $
#'
#'  Untested cases in ppp() or associated code

require(spatstat)
local({
  X <- runifpoint(10, letterR)
  Y <- runifpoint(3, complement.owin(letterR))

  #' test handling of points out-of-bounds
  df <- rbind(as.data.frame(X), as.data.frame(Y))
  A <- ppp(df$x, df$y, window=letterR, marks=1:13)
  #' test handling of points with bad coordinates
  B <- ppp(X$x, c(X$y[1:7], c(Inf, NA, NaN)), window=letterR, marks=1:10)
  D <- ppp(X$x, c(X$y[1:7], c(Inf, NA, NaN)), window=letterR,
           marks=data.frame(id=1:10, u=runif(10)))

  #' test print/summary/plot methods on these bad objects
  print(A)
  print(B)
  print(D)
  summary(A)
  summary(B)
  summary(D)
  plot(A)
  plot(B)
  plot(D)
  
  #' subset operator --- cases not covered elsewhere
  #'   subset index is a logical image
  Z <- distmap(letterR, invert=TRUE)
  V <- (Z > 0.2)
  XV <- X[V]
  #'   multiple columns of marks
  fun3 <- finpines[1:3]
  #'   multiple columns of marks, one of which is a factor
  U <- finpines
  marks(U)[,2] <- factor(c(rep("A", 60), rep("B", npoints(U)-60)))
  UU <- U[1:3, drop=TRUE]

  #' test as.ppp for spatial package if it is not installed
  FR <- Frame(letterR)
  as.ppp(list(x=X$x, y=X$y,
              xl=FR$xrange[1], xu=FR$xrange[2],
              yl=FR$yrange[1], yu=FR$yrange[2]))

  #' various utilities
  periodify(cells, 2)
  periodify(demopat, 2)

  #'
  a <- multiplicity(finpines)
  a <- multiplicity(longleaf)

  ## superimpose.ppp, extra cases
  X <- runifpoint(20)
  A <- superimpose(cells, X, W="convex")
  A <- superimpose(cells, X, W=ripras)
  ## superimpose.splitppp
  Y <- superimpose(split(amacrine))
})
#
# tests/ppx.R
#
# Test operations for ppx objects
#
#  $Revision: 1.5 $ $Date: 2019/01/02 07:58:20 $
#

require(spatstat)

local({
  #' general tests
  df <- data.frame(x=c(1,2,2,1)/4, y=c(1,2,3,1)/4, z=c(2,3,4,3)/5)
  X <- ppx(data=df, coord.type=rep("s", 3), domain=box3())
  unique(X)
  duplicated(X)
  anyDuplicated(X)
  multiplicity(X)
  print(X)
  summary(X)
  plot(X)
  domain(X)
  unitname(X) <- c("metre", "metres")
  unitname(X)

  #' subset operator
  X[integer(0)]
  Y <- X %mark% data.frame(a=df$x, b=1:4)
  Y[1:2]
  Y[FALSE]
  marks(Y) <- as.data.frame(marks(Y))
  Y[integer(0)]
  Y[1:2]
  Y[FALSE]

  #' two dimensional
  A <- ppx(data=df[,1:2], coord.type=rep("s", 2), domain=square(1))
  plot(A)
  B <- ppx(data=df[,1:2], coord.type=rep("s", 2), domain=NULL)
  plot(B)
  #' one dimensional
  E <- ppx(data=data.frame(x=runif(10)))
  plot(E)
  
  #' bug
  stopifnot(identical(unmark(chicago[1]),
                      unmark(chicago)[1]))

  #' ppx with zero points
  U <- chicago[integer(0)]
  V <- U %mark% 1
  V <- U %mark% factor("a")

  #' simplify lower-dimensional patterns
  X3 <- ppx(data=df, coord.type=rep("s", 3), domain=box3(), simplify=TRUE)
  stopifnot(is.pp3(X3))
  X2 <- ppx(data=df[,1:2], coord.type=rep("s", 2), domain=square(1), simplify=TRUE)
  stopifnot(is.ppp(X2))

  #' marks<-.ppx
  M <- as.matrix(X)
  marks(X) <- df[,1]
  marks(X) <- df[,integer(0)]

  #' trivial cases of random generators
  B4 <- boxx(0:1, 0:1, 0:1, 0:1)
  Z0 <- runifpointx(0, domain=B4, nsim=2)
  Z1 <- runifpointx(1, domain=B4, nsim=2)
})
#
# tests/prediction.R
#
# Things that might go wrong with predict()
#
#  $Revision: 1.11 $ $Date: 2019/04/27 09:25:59 $
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
  #' adaptcoef     
  a <- c(A=1,B=2,Z=42)
  b <- c(B=41,A=0)
  ab <- adaptcoef(a, b, drop=TRUE)
  
  # tests of relrisk.ppm
  fut <- ppm(amacrine ~ x * marks)
  a <- relrisk(fut, control=2, relative=TRUE)
  a <- relrisk(fut, se=TRUE)
  a <- relrisk(fut, relative=TRUE, se=TRUE)
  fut <- ppm(sporophores ~ marks + x)
  a <- relrisk(fut, control=2, relative=TRUE)
  a <- relrisk(fut, se=TRUE)
  a <- relrisk(fut, relative=TRUE, se=TRUE)

  ## untested cases of predict.ppm
  fit0 <- ppm(cells)
  a <- predict(fit0, interval="confidence")
  a <- predict(fit0, interval="confidence", type="count")
  fit  <- ppm(cells ~ x)
  b <- predict(fit, type="count",                            se=TRUE)
  b <- predict(fit, type="count", window=square(0.5),        se=TRUE)
  b <- predict(fit, type="count", window=quadrats(cells, 3), se=TRUE)
  ## superseded usages
  b <- predict(fit, type="se", getoutofjail=TRUE)
  b <- predict(fit, total=TRUE)
  b <- predict(fit, total=square(0.5))
  b <- predict(fit, total=quadrats(cells, 3))

  ## supporting code
  u <- model.se.image(fit, square(0.5))
  u <- model.se.image(fit, square(0.5), what="cv")
  u <- model.se.image(fit, square(0.5), what="ce")

  ##
  fut <- ppm(cells ~ x, Strauss(0.1))
  df <- data.frame(x=runif(10), y=runif(10),
                   Interaction=sample(0:1, 10, TRUE))
  m10 <- PPMmodelmatrix(fut, data=df)
  mmm <- PPMmodelmatrix(fut, Q=quad.ppm(fut))
  #' effectfun for Gibbs
  effectfun(fut, "x")
  effectfun(fut, "x", se.fit=TRUE)
  #' implicit covariate when there is only one
  effectfun(fut)
  effectfun(fut, se.fit=TRUE)
  #' 
  dlin <- distfun(copper$SouthLines)
  copfit <- ppm(copper$SouthPoints ~ dlin, Geyer(1,1))
  effectfun(copfit, "dlin")
  effectfun(copfit)
  # external covariate
  effectfun(fut, "y", x=0)
  futS <- ppm(cells ~ 1, Strauss(0.1))
  effectfun(futS, "x")
  effectfun(futS, "y")

  ## ppm with covariate values in data frame
  X <- rpoispp(42)
  Q <- quadscheme(X)
  weirdfunction <- function(x,y){ 10 * x^2 + 5 * sin(10 * y) }
  Zvalues <- weirdfunction(x.quad(Q), y.quad(Q))
  fot <- ppm(Q ~ y + Z, data=data.frame(Z=Zvalues))
  effectfun(fot, "y", Z=0)
  effectfun(fot, "Z", y=0)
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

  #' code coverage
  op <- spatstat.options(project.fast=TRUE)
  fut <- emend(fit, trace=TRUE)
  chk(fut)
  spatstat.options(op)

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

reset.spatstat.options()
