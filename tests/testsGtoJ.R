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
#  $Revision: 1.7 $  $Date: 2020/01/11 04:55:17 $
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

  rn <- rownames(h)
  r.n <- row.names(h)
  if(!identical(rn, r.n))
    stop("rownames and row.names conflict for hyperframes")

  dn <- dimnames(h)
  dimnames(h) <- dn
  dimnames(h)[[2]][2] <- "copacetic"
  dimnames(h)[[1]][2] <- "second"

  #' hyperframe with a hyperatom
  H <- hyperframe(A=runif(3), B=1:3, D=runifpoint(10))
  H[,3]
  H[,3,drop=TRUE]
  #' special cases of [<-
  H$B <- H[,1]
  H[2:3,1] <- H[2:3,2]
  H[2:3,1] <- H[2,2]
  H[2,1:2] <- H[3,1:2]

  #' split
  f <- factor(c("a", "a", "b"))
  G <- split(H, f)
  G[["a"]]$B <- 42
  split(H, f) <- G
})


#'     tests/hypotests.R
#'     Hypothesis tests
#' 
#'  $Revision: 1.5 $ $Date: 2019/12/14 03:11:05 $

require(spatstat)
local({
  hopskel.test(redwood, method="MonteCarlo", nsim=5)
  
  berman.test(spiders, "x")
  berman.test(lppm(spiders ~ x), "y")

  #' quadrat test - spatial methods
  a <- quadrat.test(redwood, 3)
  domain(a)
  shift(a, c(1,1))

  #' cases of studpermu.test
  #' X is a hyperframe
  b <- studpermu.test(pyramidal, nperm=9)
  b <- studpermu.test(pyramidal, nperm=9, use.Tbar=TRUE)
  #' X is a list of lists of ppp
  ZZ <- split(pyramidal$Neurons, pyramidal$group)
  bb <- studpermu.test(ZZ, nperm=9)

  #' Issue #115
  X <- runifpoint(50, nsim = 3)
  Y <- runifpoint(3000, nsim = 3)
  h <- hyperframe(ppp = c(X, Y), group = rep(1:2, 3))
  studpermu.test(h, ppp ~ group)

  #' scan test
  Z <- scanmeasure(cells, 0.1, method="fft")
  rr <- c(0.05, 1)
  scan.test(amacrine, rr, nsim=5,
            method="binomial", alternative="less")
  fit <- ppm(cells ~ x)
  lam <- predict(fit)
  scan.test(cells, rr, nsim=5,
            method="poisson", baseline=fit, alternative="less")
  scan.test(cells, rr, nsim=5,
            method="poisson", baseline=lam, alternative="less")
})
#
#  tests/imageops.R
#
#   $Revision: 1.23 $   $Date: 2020/01/21 04:52:56 $
#

require(spatstat)
local({
  #' cases of 'im' data
  tab <- table(sample(factor(letters[1:10]), 30, replace=TRUE))
  b <- im(tab, xrange=c(0,1), yrange=c(0,10))
  b <- update(b)

  mat <- matrix(sample(0:4, 12, replace=TRUE), 3, 4)
  levels(mat) <- 0:4
  b <- im(mat)
  b <- update(b)

  D <- as.im(mat, letterR)
  df <- as.data.frame(D)
  DD <- as.im(df, step=c(D$xstep, D$ystep))
  
  #' various manipulations
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
  d[,] <- d[] + 1
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

  #' cases of [.im
  Ma <- as.mask(M, dimyx=37)
  ZM <- Z[raster=Ma, drop=FALSE]
  ZM[solutionset(Y+Z > 0.4)] <- NA
  ZF <- cut(ZM, breaks=5)
  ZL <- (ZM > 0)
  P <- list(x=c(0.511, 0.774, 0.633, 0.248, 0.798),
            y=c(0.791, 0.608, 0.337, 0.613, 0.819))
  zmp <- ZM[P, drop=TRUE]
  zfp <- ZF[P, drop=TRUE]
  zlp <- ZL[P, drop=TRUE]
  P <- as.ppp(P, owin())
  zmp <- ZM[P, drop=TRUE]
  zfp <- ZF[P, drop=TRUE]
  zlp <- ZL[P, drop=TRUE]

  #' miscellaneous
  ZZ <- zapsmall(Z, digits=6)
  ZZ <- zapsmall(Z)

  ZS <- shift(Z, origin="centroid")
  ZS <- shift(Z, origin="bottomleft")

  ZA <- affine(Z, mat=diag(c(-1,-2)))
  
  #' hist.im
  h <- hist(Z)
  h <- hist(Z, plot=FALSE)
  Zcut <- cut(Z, breaks=5)
  h <- hist(Zcut) # barplot
  plot(h) # plot.barplotdata

  #' plot.im code blocks
  plot(Z, ribside="left")
  plot(Z, ribside="top")
  plot(Z, riblab="value")
  plot(Z, clipwin=square(0.5))
  plot(Z - mean(Z), log=TRUE)
  plot(Z, valuesAreColours=TRUE) # rejected with a warning
  IX <- as.im(function(x,y) { as.integer(round(3*x)) }, square(1))
  co <- colourmap(rainbow(4), inputs=0:3)
  plot(IX, col=co)
  CX <- eval.im(col2hex(IX+1L))
  plot(CX, valuesAreColours=TRUE)
  plot(CX, valuesAreColours=FALSE)

  #' pairs.im 
  pairs(solist(Z))
  pairs(solist(A=Z))
  
  #' handling and plotting of character and factor images
  Afactor    <- as.im(col2hex("green"), letterR, na.replace=col2hex("blue"))
  Acharacter <- as.im(col2hex("green"), letterR, na.replace=col2hex("blue"),
                      stringsAsFactors=FALSE)
  plot(Afactor)
  plot(Acharacter, valuesAreColours=TRUE)
  print(summary(Afactor))
  print(summary(Acharacter))

  #' safelookup (including extrapolation case)
  Z <- as.im(function(x,y) { x - y }, letterR)
  B <- grow.rectangle(Frame(letterR), 1)
  X <- superimpose(runifpoint(10,letterR),
                   runifpoint(20, setminus.owin(B, letterR)),
                   W=B)
  a <- safelookup(Z, X)

  #' check nearest.valid.pixel
  W <- Window(demopat)
  set.seed(911911)
  X <- runifpoint(1000, W)
  Z <- quantess(W, function(x,y) { x }, 9)$image
  x <- X$x
  y <- X$y
  a <- nearest.valid.pixel(x, y, Z, method="interpreted")
  b <- nearest.valid.pixel(x, y, Z, method="C")
  if(!isTRUE(all.equal(a,b)))
    stop("Unequal results in nearest.valid.pixel")
  if(!identical(a,b)) 
    stop("Equal, but not identical, results in nearest.valid.pixel")

  #' cases of distcdf
  distcdf(cells[1:5])
  distcdf(W=cells[1:5], dW=1:5)
  distcdf(W=Window(cells), V=cells[1:5])
  distcdf(W=Window(cells), V=cells[1:5], dV=1:5)

  #' im.apply
  DA <- density(split(amacrine))
  Z <- im.apply(DA, sd)
  Z <- which.max.im(DA) # deprecated -> im.apply(DA, which.max)

  #' Math.imlist, Ops.imlist, Complex.imlist
  U <- Z+2i
  B <- U * (2+1i)
  print(summary(B))
  V <- solist(A=U, B=B)
  negV <- -V
  E <- Re(V)
  negE <- -E

  #' rotmean
  U <- rotmean(Z, origin="midpoint", result="im", padzero=FALSE)
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
#'
#'  tests/interact.R
#'
#'  Support for interaction objects
#'
#'  $Revision: 1.1 $ $Date: 2019/12/10 01:57:18 $

require(spatstat)
local({
  #' print.intermaker
  Strauss
  Geyer
  Ord
  #' intermaker
  BS <- get("BlankStrauss", envir=environment(Strauss))
  BD <- function(r) { instantiate.interact(BS, list(r=r)) }
  BlueDanube <- intermaker(BD, BS) 
})
#'   tests/ippm.R
#'   Tests of 'ippm' class
#'   $Revision: 1.4 $ $Date: 2020/01/07 09:36:42 $

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
  # fit model with logistic likelihood but without iScore
  fitlo <- ippm(X ~Z + offset(log(f)),
                method="logi",
                covariates=list(Z=Z, f=f),
                start=list(gamma=1, delta=1),
                nd=nd)

  ## ............. test ippm class support ......................
  Ar <- model.matrix(fit)
  Ai <- model.matrix(fit, irregular=TRUE)
  Zr <- model.images(fit)
  Zi <- model.images(fit, irregular=TRUE)
  ## update.ippm
  fit2 <- update(fit, . ~ . + I(Z^2))
  fit0 <- update(fit,
                 . ~ . - Z,
                 start=list(gamma=2, delta=4))
  oldfit <- ippm(X,
              ~Z + offset(log(f)),
              covariates=list(Z=Z, f=f),
              iScore=Dlogf,
              start=list(gamma=1, delta=1),
              nd=nd)
  oldfit2 <- update(oldfit, . ~ . + I(Z^2))
  oldfit0 <- update(oldfit,
                    . ~ . - Z,
                    start=list(gamma=2, delta=4))
  ## again with logistic
  fitlo2 <- update(fitlo, . ~ . + I(Z^2))
  fitlo0 <- update(fitlo,
                   . ~ . - Z,
                   start=list(gamma=2, delta=4))
  oldfitlo <- ippm(X,
                   ~Z + offset(log(f)),
                   method="logi",
                   covariates=list(Z=Z, f=f),
                   start=list(gamma=1, delta=1),
                   nd=nd)
  oldfitlo2 <- update(oldfitlo, . ~ . + I(Z^2))
  oldfitlo0 <- update(oldfitlo,
                      . ~ . - Z,
                      start=list(gamma=2, delta=4))
  ## anova.ppm including ippm objects
  fit0 <- update(fit, . ~ Z)
  fit0lo <- update(fitlo, . ~ Z)
  A <- anova(fit0, fit)
  Alo <- anova(fit0lo, fitlo)
})
