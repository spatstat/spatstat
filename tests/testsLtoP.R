## 
##    tests/legacy.R
##
## Test that current version of spatstat is compatible with outmoded usage
## $Revision: 1.2 $ $Date: 2015/12/29 08:54:49 $
local({
  require(spatstat)

  ## (1) Old syntax of ppm
  ppm(cells, ~x)
  
  ## (2) Old syntax of MultiStrauss etc.
  r <- matrix(3, 2, 2)
  a <- MultiStrauss( , r)
  a <- MultiStrauss(NULL, r)
  a <- MultiHard(, r)
  
  h <- r/2
  a <- MultiStraussHard( , r, h)

  NULL
})
#'
#'    tests/leverinf.R
#'
#'   leverage and influence for Gibbs models
#' 
#'   $Revision: 1.5 $ $Date: 2016/07/06 08:29:49 $
#' 

require(spatstat)
local({
  # original non-sparse algorithm
  Leverage <- function(...) leverage(..., sparseOK=FALSE)
  Influence <- function(...) influence(..., sparseOK=FALSE)
  Dfbetas <- function(...) dfbetas(..., sparseOK=FALSE)
  # Strauss()$delta2
  fitS <- ppm(cells ~ x, Strauss(0.12), rbord=0)
  levS <- Leverage(fitS)
  infS <- Influence(fitS)
  dfbS <- Dfbetas(fitS)
  # Geyer()$delta2
  fitG <- ppm(redwood ~ 1, Geyer(0.1, 2), rbord=0)
  levG <- Leverage(fitG)
  infG <- Influence(fitG)
  # pairwise.family$delta2
  fitD <- ppm(cells ~ 1, DiggleGatesStibbard(0.12), rbord=0)
  levD <- Leverage(fitD)
  infD <- Influence(fitD)
  # ppmInfluence; offset is present; coefficient vector has length 1
  fitH <- ppm(cells ~ x, Hardcore(0.07), rbord=0)
  levH <- Leverage(fitH)
  infH <- Influence(fitH)

  ## divide and recombine algorithm
  op <- spatstat.options(maxmatrix=50000)
  levSB <- Leverage(fitS)
  infSB <- Influence(fitS)
  dfbSB <- Dfbetas(fitS)

  chk <- function(x, y, what,
                  from="single-block and multi-block",
                  thresh=1e-12) {
    if(max(abs(x-y)) > thresh)
      stop(paste("Different results for", what, "obtained from",
                 from, "algorithms"),
           call.=FALSE)
    invisible(NULL)
  }

  chk(marks(as.ppp(infS)), marks(as.ppp(infSB)), "influence")
  chk(as.im(levS),         as.im(levSB),         "leverage")
  chk(dfbS$val,            dfbSB$val,            "dfbetas$value")
  chk(dfbS$density,        dfbSB$density,        "dfbetas$density")

  spatstat.options(op)

  ## sparse algorithm
  pmi <- ppmInfluence(fitS, sparseOK=TRUE)
  levSp <- pmi$leverage
  infSp <- pmi$influence
  dfbSp <- pmi$dfbetas

  chks <- function(...) chk(..., from="sparse and non-sparse")
  
  chks(marks(as.ppp(infS)), marks(as.ppp(infSp)), "influence")
  chks(as.im(levS),         as.im(levSp),         "leverage")
  chks(dfbS$val,            dfbSp$val,            "dfbetas$value")
  chks(dfbS$density,        dfbSp$density,        "dfbetas$density")
})

##
##    tests/linalgeb.R
##
## checks validity of linear algebra code
##
##  $Revision: 1.3 $ $Date: 2015/12/29 08:54:49 $
##
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
## 
## tests/localpcf.R
##
## temporary test file for localpcfmatrix
##  $Revision: 1.2 $  $Date: 2015/12/29 08:54:49 $

require(spatstat)
local({
  a <- localpcfmatrix(redwood)
  a
  plot(a)
  a[, 3:5]
})
#
# tests/lppstuff.R
#
# Tests for lpp code
#
#  $Revision: 1.8 $  $Date: 2016/07/18 06:43:57 $


require(spatstat)

local({
  # check 'normalise' option in linearKinhom
  X <- rpoislpp(5, simplenet)
  fit <- lppm(X ~x)
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
  
  ## nearest neighbour distances
  eps <- sqrt(.Machine$double.eps)
  f <- function(mat,k) { apply(mat, 1, function(z,n) { sort(z)[n]  }, n=k+1) }
  g <- function(mat,k) { apply(mat, 1, function(z,n) { order(z)[n] }, n=k+1) }

  XX <- spiders
  nn <- nndist(XX)
  nnP <- f(pairdist(XX), 1)
  if(any(abs(nn - nnP) > eps))
    stop("nndist.lpp does not agree with pairdist.lpp")

  nw <- nnwhich(XX)
  nwP <- g(pairdist(XX), 1)
  if(any(nw != nwP))
    stop("nnwhich.lpp does not agree with pairdist")

  ZZ <- split(chicago)
  XX <- ZZ$damage
  YY <- ZZ$assault
  op <- spatstat.options(Cnncrosslpp=FALSE)
  a <- nncross(XX, YY)
  spatstat.options(Cnncrosslpp=TRUE)
  b <- nncross(XX, YY)
  if(any(a$which != b$which))
    stop("Inconsistent values of nncross.lpp()$which from different C code")
  if(max(abs(a$dist - b$dist)) > eps)
    stop("Inconsistent values of nncross.lpp()$dist from different C code")

  spatstat.options(Cnncrosslpp=TRUE)
  b2 <- nncross(XX, YY, k=1:2, what="which")
  if(any(b2$which.1 != b$which))
    stop("inconsistent values of nncross.lpp()$which from k=1:2 and k=1")
  a2 <- nncross(XX, YY, k=1:2, what="dist")
  if(max(abs(a2$dist.1 - a$dist)) > eps)
    stop("Inconsistent values of nncross.lpp()$dist from k=1:2 and k=1")

  spatstat.options(Cnncrosslpp=TRUE)
  ii <- seq_len(npoints(XX))
  w1 <- nnwhich(XX)
  w2 <- nncross(XX, XX, iX=ii, iY=ii, what="which")
  w3 <- nncross(XX, XX, iX=ii, iY=ii, what="which", method="interpreted")
  if(any(w1 != w2))
    stop("nnwhich.lpp disagrees with nncross.lpp(iX, iY)")
  if(any(w2 != w3))
    stop("Different results for nncross.lpp(iX, iY, 'which') using R and C")
  d1 <- nndist(XX)
  d2 <- nncross(XX, XX, iX=ii, iY=ii, what="dist")
  d3 <- nncross(XX, XX, iX=ii, iY=ii, what="dist", method="interpreted")
  if(max(abs(d1-d2)) > eps)
    stop("nndist.lpp disagrees with nncross.lpp(iX, iY)")
  if(max(abs(d2-d3)) > eps)
    stop("Different results for nncross.lpp(iX, iY, 'dist') using R and C")

  spatstat.options(op)

  # test handling marginal cases
  xyd <- nncross(XX, YY[1])

  ## Test linnet surgery code
  set.seed(42)
  X <- runiflpp(30, simplenet)
  V <- runiflpp(30, simplenet)
  XV <- insertVertices(X, V)
  validate.lpp.coords(XV, context="calculated by insertVertices")

  ## Test [.lpp internal data
  B <- owin(c(0.1,0.7),c(0.19,0.6))
  XB <- X[B]
  validate.lpp.coords(XB, context="returned by [.lpp")

  ## Tests related to linearK, etc
  testcountends <- function(X, r=100, s=1) {
    if(s != 1) {
      X <- rescale(X, s)
      r <- r/s
    }
    L <- as.linnet(X)
    n1 <- countends(L, X[1], r)
    n2 <- npoints(lineardisc(L, X[1], r, plotit=FALSE)$endpoints)
    if(n1 != n2)
      stop(paste("Incorrect result from countends:",
                 n1, "!=", n2, 
                 paren(paste("scale=", 1/s))),
           call.=FALSE)
  }
  # original scale
  X <- unmark(chicago)
  testcountends(X)
  # finer scale
  testcountends(X, s=1000)

  ## Test algorithms for boundingradius.linnet
  L <- as.linnet(chicago, sparse=TRUE)
  opa <- spatstat.options(Clinearradius=FALSE)
  bR <- as.linnet(L, sparse=FALSE)$boundingradius
  spatstat.options(Clinearradius=TRUE)
  bC <- as.linnet(L, sparse=FALSE)$boundingradius
  spatstat.options(opa)
  if(abs(bR-bC) > 0.001 * (bR+bC)/2)
    stop("Disagreement between R and C algorithms for boundingradius.linnet",
         call.=FALSE)
})

##
##     tests/marcelino.R
##
##     $Revision: 1.3 $  $Date: 2015/12/29 08:54:49 $
##
require(spatstat)

local({
  Y <- split(urkiola)
  B <- Y$birch
  O <- Y$oak
  B.lam <- predict (ppm(B ~polynom(x,y,2)), type="trend")
  O.lam <- predict (ppm(O ~polynom(x,y,2)), type="trend")

  Kinhom(B, lambda=B.lam, correction="iso")
  Kinhom(B, lambda=B.lam, correction="border")

  Kcross.inhom(urkiola, i="birch", j="oak", B.lam, O.lam)
  Kcross.inhom(urkiola, i="birch", j="oak", B.lam, O.lam, correction = "iso")
  Kcross.inhom(urkiola, i="birch", j="oak", B.lam, O.lam, correction = "border")
})


##
##    tests/markcor.R
##
##   Tests of mark correlation code (etc)
##
## $Revision: 1.4 $ $Date: 2015/12/29 08:54:49 $

require(spatstat)

local({
  ## check.testfun checks equality of functions
  ##  and is liable to break if the behaviour of all.equal is changed
  fe <- function(m1, m2) {m1 == m2}
  fm <- function(m1, m2) {m1 * m2}
  fs <- function(m1, m2) {sqrt(m1)}
  if(check.testfun(fe, X=amacrine)$ftype != "equ")
    warning("check.testfun fails to recognise mark equality function")
  if(check.testfun(fm, X=longleaf)$ftype != "mul")
    warning("check.testfun fails to recognise mark product function")
  check.testfun(fs, X=longleaf)
  
  ## test all is well in Kmark -> Kinhom 
  MA <- Kmark(amacrine,function(m1,m2){m1==m2})
  set.seed(42)
  AR <- rlabel(amacrine)
  MR <- Kmark(AR,function(m1,m2){m1==m2})
  if(isTRUE(all.equal(MA,MR)))
    stop("Kmark unexpectedly ignores marks")
})
#
# tests/mppm.R
#
# Basic tests of mppm
#
# $Revision: 1.8 $ $Date: 2016/06/28 04:19:08 $
# 

require(spatstat)

local({
  ## test interaction formulae and subfits
  fit1 <- mppm(Points ~ group, simba,
               hyperframe(po=Poisson(), str=Strauss(0.1)),
               iformula=~ifelse(group=="control", po, str))
  fit2 <- mppm(Points ~ group, simba,
               hyperframe(po=Poisson(), str=Strauss(0.1)),
               iformula=~str/id)
  fit3 <- mppm(Points ~ group, simba,
               hyperframe(po=Poisson(), pie=PairPiece(c(0.05,0.1))),
        iformula=~I((group=="control") * po) + I((group=="treatment") * pie))
  fit1
  fit2
  fit3

  ## run summary.mppm which currently sits in spatstat-internal.Rd
  summary(fit1)
  summary(fit2)
  summary(fit3)

  ## test vcov algorithm
  vcov(fit1)
  vcov(fit2)
  vcov(fit3)

  ## test subfits algorithm
  s1 <- subfits(fit1)
  s2 <- subfits(fit2)
  s3 <- subfits(fit3)

  ## validity of results of subfits()
  p1 <- solapply(s1, predict)
  p2 <- solapply(s2, predict)
  p3 <- solapply(s3, predict)

})

local({
  ##  [thanks to Sven Wagner]
  ## factor covariate, with some levels unused in some rows
  set.seed(14921788)
  H <- hyperframe(X=replicate(3, runifpoint(20), simplify=FALSE),
                  Z=solist(as.im(function(x,y){x}, owin()),
                    as.im(function(x,y){y}, owin()),
                    as.im(function(x,y){x+y}, owin())))
  H$Z <- solapply(H$Z, cut, breaks=(0:4)/2)

  fit6 <- mppm(X ~ Z, H)
  v6 <- vcov(fit6)
  s6 <- subfits(fit6)
  p6 <- solapply(s6, predict)

  # random effects
  fit7 <- mppm(X ~ Z, H, random=~1|id)
  v7 <- vcov(fit7)
  s7 <- subfits(fit7)
  p7 <- solapply(s7, predict)

  fit7a <- mppm(X ~ Z, H, random=~x|id)
  v7a <- vcov(fit7a)
  s7a <- subfits(fit7a)
  p7a <- solapply(s7a, predict)

  # multitype: collisions in vcov.ppm, predict.ppm
  H$X <- lapply(H$X, rlabel, labels=factor(c("a","b")), permute=FALSE)
  M <- MultiStrauss(matrix(0.1, 2, 2), c("a","b"))
  fit8 <- mppm(X ~ Z, H, M)
  v8 <- vcov(fit8, fine=TRUE)
  s8 <- subfits(fit8)
  p8 <- lapply(s8, predict)
  c8 <- lapply(s8, predict, type="cif")

  fit9 <- mppm(X ~ Z, H, M, iformula=~Interaction * id)
  v9 <- vcov(fit9, fine=TRUE)
  s9 <- subfits(fit9)
  p9 <- lapply(s9, predict)
  c9 <- lapply(s9, predict, type="cif")

  # and a simple error in recognising 'marks'
  fit10 <- mppm(X ~ marks, H)
})

local({
  ## test handling of offsets and zero cif values in mppm
  H <- hyperframe(Y = waterstriders)
  mppm(Y ~ 1,  data=H, Hardcore(1.5))
  mppm(Y ~ 1,  data=H, StraussHard(7, 1.5))

  ## prediction, in training/testing context
  ##    (example from Markus Herrmann and Ege Rubak)
  X <- waterstriders
  dist <- solapply(waterstriders,
                   function(z) distfun(runifpoint(1, Window(z))))
  i <- 3
  train <- hyperframe(pattern = X[-i], dist = dist[-i])
  test <- hyperframe(pattern = X[i], dist = dist[i])
  fit <- mppm(pattern ~ dist, data = train)
  pred <- predict(fit, type="cif", newdata=test, verbose=TRUE)
})

local({
  ## test handling of interaction coefficients in multitype case
  set.seed(42)
  XX <- as.solist(replicate(3, rthin(amacrine, 0.8), simplify=FALSE))
  H <- hyperframe(X=XX)
  M <- MultiStrauss(matrix(0.1, 2, 2), levels(marks(amacrine)))
  fit <- mppm(X ~ 1, H, M)
  co <- coef(fit)
  subco <- sapply(subfits(fit), coef)
  if(max(abs(subco - co)) > 0.001)
    stop("Wrong coefficient values in subfits, for multitype interaction")
})
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
#   $Revision: 1.16 $  $Date: 2015/12/29 08:54:49 $
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
