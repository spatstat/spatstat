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
#  $Revision: 1.5 $  $Date: 2019/09/11 11:33:37 $
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
})


#'     tests/hypotests.R
#'     Hypothesis tests
#' 
#'  $Revision: 1.4 $ $Date: 2019/10/16 02:38:39 $

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
})
#
#  tests/imageops.R
#
#   $Revision: 1.20 $   $Date: 2019/12/03 07:29:10 $
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

  #' handling and plotting of character and factor images
  Afactor    <- as.im(col2hex("green"), letterR, na.replace=col2hex("blue"))
  Acharacter <- as.im(col2hex("green"), letterR, na.replace=col2hex("blue"),
                      stringsAsFactors=FALSE)
  plot(Afactor)
  plot(Acharacter, valuesAreColours=TRUE)

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
#'   $Revision: 1.3 $ $Date: 2019/02/15 10:08:26 $

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

})
#'
#'   tests/Kfuns.R
#'
#'   Various K and L functions and pcf
#'
#'   $Revision: 1.30 $  $Date: 2019/12/08 02:19:15 $
#'

require(spatstat)
myfun <- function(x,y){(x+1) * y } # must be outside
local({
  #' supporting code
  implemented.for.K(c("border", "bord.modif", "translate", "good", "best"),
                    "polygonal", TRUE)
  implemented.for.K(c("border", "bord.modif", "translate", "good", "best"),
                    "mask", TRUE)
  implemented.for.K(c("border", "isotropic"), "mask", FALSE)
  #' Kest special code blocks
  K <- Kest(cells, var.approx=TRUE, ratio=FALSE)
  Z <- distmap(cells) + 1
  Kb <- Kest(cells, correction=c("border","bord.modif"),
             weights=Z, ratio=TRUE)
  Kn <- Kest(cells, correction="none",
             weights=Z, ratio=TRUE)
  Knb <- Kest(cells, correction=c("border","bord.modif","none"),
              weights=Z, ratio=TRUE)
  bigint <- 50000 # This is only "big" on a 32-bit system where
                  # sqrt(.Machine$integer.max) = 46340.9
  X <- runifpoint(bigint)
  Z <- as.im(1/bigint, owin())
  Kb <- Kest(X, correction=c("border","bord.modif"),
             rmax=0.02, weights=Z, ratio=TRUE)
  Kn <- Kest(X, correction="none",
             rmax=0.02, weights=Z, ratio=TRUE)
  Knb <- Kest(X, correction=c("border","bord.modif","none"),
              rmax=0.02, weights=Z, ratio=TRUE)
  #' pcf.ppp special code blocks
  pr  <- pcf(cells, ratio=TRUE, var.approx=TRUE)
  pc  <- pcf(cells, domain=square(0.5))
  pcr <- pcf(cells, domain=square(0.5), ratio=TRUE)
  pw <- pcf(redwood, correction="none")
  pwr <- pcf(redwood, correction="none", ratio=TRUE)
  #' Kinhom code blocks
  X <- rpoispp(function(x,y) { 100 * x }, 100, square(1))
  lambda <- 100 * X$x
  Kin <- Kinhom(X, lambda, correction=c("none", "border"))
  lambda2 <- outer(lambda, lambda, "*")
  Ki2 <- Kinhom(X, lambda2=lambda2, diagonal=FALSE,
                correction=c("translate", "isotropic"))
  fut <- ppm(X ~ x)
  Kio <- Kinhom(X, fut, update=FALSE)
  Kiu <- Kinhom(X, fut, update=TRUE, diagonal=FALSE)
  #' inhomogeneous multitype
  fit <- ppm(amacrine ~ marks)
  K1 <- Kcross.inhom(amacrine, lambdaX=fit)
  K2 <- Kcross.inhom(amacrine, lambdaX=densityfun(amacrine))
  K3 <- Kcross.inhom(amacrine, lambdaX=density(amacrine, at="points"))
  On <- split(amacrine)$on
  Off <- split(amacrine)$off
  K4 <- Kcross.inhom(amacrine, lambdaI=ppm(On), lambdaJ=ppm(Off))
  K5 <- Kcross.inhom(amacrine, correction="bord.modif")
  #' Kmark (=markcorrint)
  X <- runifpoint(100) %mark% runif(100)
  km <- Kmark(X, f=atan2)
  km <- Kmark(X, f1=sin)
  km <- Kmark(X, f="myfun")
  aa <- Kmark(X, normalise=FALSE, returnL=FALSE)
  aa <- Kmark(X, normalise=FALSE, returnL=TRUE)
  aa <- Kmark(X, normalise=TRUE,  returnL=FALSE)
  aa <- Kmark(X, normalise=TRUE,  returnL=TRUE)
  #'
  rr <- rep(0.1, npoints(cells))
  eC <- edge.Ripley(cells, rr)
  eI <- edge.Ripley(cells, rr, method="interpreted")
  if(max(abs(eC-eI)) > 0.1)
    stop("Ripley edge correction results do not match")

  a <- rmax.Ripley(square(1))
  a <- rmax.Rigid(square(1))
  a <- rmax.Ripley(as.polygonal(square(1)))
  a <- rmax.Rigid(as.polygonal(square(1)))
  a <- rmax.Ripley(letterR)
  a <- rmax.Rigid(letterR)

  #' run slow code for edge correction and compare results
  X <- redwood[c(TRUE, FALSE, FALSE)]
  Window(X) <- as.polygonal(Window(X))
  Eapprox <- edge.Trans(X)
  Eexact <- edge.Trans(X, exact=TRUE)
  maxrelerr <- max(abs(1 - range(Eapprox/Eexact)))
  if(maxrelerr > 0.1)
    stop(paste("Exact and approximate algorithms for edge.Trans disagree by",
               paste0(round(100*maxrelerr), "%")),
         call.=FALSE)
  #'
  #'   directional K functions
  #'
    a <- Ksector(swedishpines,
                 -pi/2, pi/2, units="radians",
                 correction=c("none", "border", "bord.modif", "Ripley", "translate"),
                 ratio=TRUE)
    plot(a)
                 
  #'
  #'   local K functions
  #'
  fut <- ppm(swedishpines ~ polynom(x,y,2))
  Z <- predict(fut)
  Lam <- fitted(fut, dataonly=TRUE)
  a <- localLinhom(swedishpines, lambda=fut)
  a <- localLinhom(swedishpines, lambda=Z)
  a <- localLinhom(swedishpines, lambda=Lam)
  a <- localLinhom(swedishpines, lambda=Z, correction="none")
  a <- localLinhom(swedishpines, lambda=Z, correction="translate")
  a <- localLcross(amacrine)
  a <- localKdot(amacrine)
  a <- localLdot(amacrine)
  a <- localKcross.inhom(amacrine)
  a <- localLcross.inhom(amacrine)
  fat <- ppm(amacrine ~ x * marks)
  Zed <- predict(fat)
  Lum <- fitted(fat, dataonly=TRUE)
  moff <- (marks(amacrine) == "off")
  a <- localLcross.inhom(amacrine, from="off", to="on", lambdaX=Zed)
  a <- localLcross.inhom(amacrine, from="off", to="on", lambdaX=Lum)
  a <- localLcross.inhom(amacrine, from="off", to="on", lambdaX=fat)
  a <- localLcross.inhom(amacrine, from="off", to="on",
                         lambdaFrom=Lum[moff], lambdaTo=Lum[!moff])
  a <- localLcross.inhom(amacrine, from="off", to="on", lambdaX=Zed,
                         correction="none")
  a <- localLcross.inhom(amacrine, from="off", to="on", lambdaX=Zed,
                         correction="translate")
  #'
  #' cases of resolve.lambda.cross
  #'
  h <- resolve.lambda.cross(amacrine, moff, !moff)
  h <- resolve.lambda.cross(amacrine, moff, !moff, lambdaX=Zed)
  h <- resolve.lambda.cross(amacrine, moff, !moff, lambdaX=Lum)
  h <- resolve.lambda.cross(amacrine, moff, !moff, lambdaX=fat)
  h <- resolve.lambda.cross(amacrine, moff, !moff, lambdaX=fat, update=FALSE)
  h <- resolve.lambda.cross(amacrine, moff, !moff,
                            lambdaI=Zed[["off"]], lambdaJ=Zed[["on"]])
  h <- resolve.lambda.cross(amacrine, moff, !moff,
                            lambdaI=Lum[moff], lambdaJ=Lum[!moff])
  h <- resolve.lambda.cross(amacrine, moff, !moff,
                            lambdaI=fat, lambdaJ=fat)
  h <- resolve.lambda.cross(amacrine, moff, !moff,
                            lambdaI=fat, lambdaJ=fat,
                            update=FALSE)
  #'
  #'   lohboot code blocks
  #'
  Ared <- lohboot(redwood, fun="Kest", block=TRUE,
                  Vcorrection=TRUE, global=FALSE, correction="none")
  Bred <- lohboot(redwood, block=TRUE, basicboot=TRUE, global=FALSE)
  Cred <- lohboot(redwood, fun=Kest, block=TRUE, global=TRUE,
                  correction="translate")
  Dred <- lohboot(redwood, Lest)
  Kred <- lohboot(redwood, Kinhom)
  Lred <- lohboot(redwood, Linhom)
  gred <- lohboot(redwood, pcfinhom, sigma=0.1)
  Zred <- predict(ppm(redwood ~ x+y))
  Lred <- lohboot(redwood, Linhom, lambda=Zred)
  #'
  X <- runifpoint(100, letterR)
  AX <- lohboot(X, block=TRUE, nx=7, ny=10)
  #'    multitype
  b <- lohboot(amacrine, Kcross)
  b <- lohboot(amacrine, Lcross)
  b <- lohboot(amacrine, Kdot)
  b <- lohboot(amacrine, Ldot)
  b <- lohboot(amacrine, Kcross.inhom)
  b <- lohboot(amacrine, Lcross.inhom)
  b <- lohboot(amacrine, Lcross.inhom, from="off", to="on", lambdaX=Zed)
  b <- lohboot(amacrine, Lcross.inhom, from="off", to="on", lambdaX=Lum)
  b <- lohboot(amacrine, Lcross.inhom, from="off", to="on", lambdaX=fat)
  b <- lohboot(amacrine, Lcross.inhom, from="off", to="on", 
                         lambdaFrom=Lum[moff], lambdaTo=Lum[!moff])
  #'
  #'  residual K functions etc
  #'
  rco <- compareFit(cells, Kcom,
                    interaction=anylist(P=Poisson(), S=Strauss(0.08)),
                    same="trans", different="tcom")
})
  
local({
  #' From Ege, in response to a stackoverflow question.
  #' The following example has two points separated by r = 1 with 1/4 of the
  #' circumference outside the 10x10 window (i.e. area 100).
  #' Thus the value of K^(r) should jump from 0 to 
  #' 100/(2\cdot 1)\cdot ((3/4)^{-1} + (3/4)^{-1}) = 100 \cdot 4/3 = 133.333.
  x <- c(4.5,5.5)
  y <- c(10,10)-sqrt(2)/2
  W <- square(10)
  X <- ppp(x, y, W)
  compere <- function(a, b, where, tol=1e-6) {
    descrip <- paste("discrepancy in isotropic edge correction", where)
    err <- as.numeric(a) - as.numeric(b)
    maxerr <- max(abs(err))
    blurb <- paste(descrip, "is", paste0(signif(maxerr, 4), ","), 
                   if(maxerr > tol) "exceeding" else "within",
                   "tolerance of", tol)
    message(blurb)
    if(maxerr > tol) {
      message(paste("Discrepancies:", paste(err, collapse=", ")))
      stop(paste("excessive", descrip), call.=FALSE)
    }
    invisible(TRUE)
  }
  ## Testing:
  eX <- edge.Ripley(X, c(1,1))
  compere(eX, c(4/3,4/3), "at interior point of rectangle")
  ## Corner case:
  Y <- X
  Y$x <- X$x-4.5+sqrt(2)/2
  eY <- edge.Ripley(Y, c(1,1))
  compere(eY, c(2,4/3), "near corner of rectangle")
  ## Invoke polygonal code
  Z <- rotate(Y, pi/4)
  eZdebug <- edge.Ripley(Z, c(1,1), internal=list(debug=TRUE))
  compere(eZdebug, c(2,4/3), "at interior point of polygon (debug on)")
  ## test validity without debugger, in case of quirks of compiler optimisation
  eZ <- edge.Ripley(Z, c(1,1))
  compere(eZ,      c(2,4/3), "at interior point of polygon (debug off)")
})
#
# tests/kppm.R
#
# $Revision: 1.31 $ $Date: 2019/11/28 21:43:33 $
#
# Test functionality of kppm that depends on RandomFields
# Test update.kppm for old style kppm objects

require(spatstat)
local({

 fit <- kppm(redwood ~1, "Thomas") # sic
 fitx <- kppm(redwood ~x, "Thomas", verbose=TRUE)
 fitx <- update(fit, ~ . + x)
 fitM <- update(fit, clusters="MatClust")
 fitC <- update(fit, cells)
 fitCx <- update(fit, cells ~ x)

 #'
 Wsub <- owin(c(0, 0.5), c(-0.5, 0))
 Zsub <- (bdist.pixels(Window(redwood)) > 0.1)
 fitWsub <- kppm(redwood ~1, "Thomas", subset=Wsub)
 fitZsub <- kppm(redwood ~1, "Thomas", subset=Zsub)
 fitWsub
 
 #' various methods
 ff <- as.fv(fitx)
 Y <- simulate(fitx, seed=42, saveLambda=TRUE)[[1]]
 uu <- unitname(fitx)
 unitname(fitCx) <- "furlong"
 mo <- model.images(fitCx)
 p <- psib(fit)
 px <- psib(fitx)
 
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
 fitxIs <- update(fitx, improve.type="quasi", fast=FALSE) 
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
   Y0 <- simulate(fit0, saveLambda=TRUE)[[1]]
   stopifnot(is.ppp(Y0))
   p0 <- psib(fit0) # issues a warning

   ## fit LGCP using K function: slow
   fit1 <- kppm(redwood ~x, "LGCP",
                covmodel=list(model="matern", nu=0.3),
                control=list(maxit=3))
   Y1 <- simulate(fit1, saveLambda=TRUE)[[1]]
   stopifnot(is.ppp(Y1))

   ## fit LGCP using pcf
   fit1p <- kppm(redwood ~x, "LGCP",
                 covmodel=list(model="matern", nu=0.3),
                 statistic="pcf")
   Y1p <- simulate(fit1p, saveLambda=TRUE)[[1]]
   stopifnot(is.ppp(Y1p))

   ## .. and using different fitting methods
   fit1pClik <- update(fit1p, method="clik")
   fit1pPalm <- update(fit1p, method="palm")
   
   ## image covariate (a different code block) 
   xx <- as.im(function(x,y) x, Window(redwood))
   fit1xx <- update(fit1p, . ~ xx, data=solist(xx=xx))
   Y1xx <- simulate(fit1xx, saveLambda=TRUE)[[1]]
   stopifnot(is.ppp(Y1xx))
   fit1xxVG <- update(fit1xx, clusters="VarGamma", nu=-1/4)
   Y1xxVG <- simulate(fit1xxVG, saveLambda=TRUE)[[1]]
   stopifnot(is.ppp(Y1xxVG))
   fit1xxLG <- update(fit1xx, clusters="LGCP",
                      covmodel=list(model="matern", nu=0.3),
                      statistic="pcf")
   Y1xxLG <- simulate(fit1xxLG, saveLambda=TRUE, drop=TRUE)
   stopifnot(is.ppp(Y1xxLG))
   
   # ... and Abdollah's code

   fit2 <- kppm(redwood ~x, cluster="Cauchy", statistic="K")
   Y2 <- simulate(fit2, saveLambda=TRUE)[[1]]
   stopifnot(is.ppp(Y2))

 }

})

local({
  #' various code blocks
  fut <- kppm(redwood, ~x)
  fet <- update(fut, redwood3)
  fot <- update(fut, trend=~y)
  fit <- kppm(redwoodfull ~ x)
  Y <- simulate(fit, window=redwoodfull.extra$regionII, saveLambda=TRUE)
  gut <- improve.kppm(fit, type="wclik1")
  gut <- improve.kppm(fit, vcov=TRUE, fast.vcov=TRUE, save.internals=TRUE)
  hut <- kppm(redwood ~ x, method="clik", weightfun=NULL)
  hut <- kppm(redwood ~ x, method="palm", weightfun=NULL)
  mut <- kppm(redwood)
  nut <- update(mut, Y)
})

local({
  #' minimum contrast code
  K <- Kest(redwood)
  a <- matclust.estK(K)
  a <- thomas.estK(K)
  a <- cauchy.estK(K)
  a <- vargamma.estK(K)
  a <- lgcp.estK(K)

  print(a)
  u <- unitname(a)
  
  g <- pcf(redwood)
  a <- matclust.estpcf(g)
  a <- thomas.estpcf(g)
  a <- cauchy.estpcf(g)
  a <- vargamma.estpcf(g)
  a <- lgcp.estpcf(g)

  #' auxiliary functions
  b <- resolve.vargamma.shape(nu.pcf=1.5)
  Z <- clusterfield("Thomas", kappa=1, scale=0.2)
  
  aa <- NULL
  aa <- accumulateStatus(simpleMessage("Woof"), aa)
  aa <- accumulateStatus(simpleMessage("Sit"),  aa)
  aa <- accumulateStatus(simpleMessage("Woof"), aa)
  printStatusList(aa)

  RMIN <- 0.01
  fit <- kppm(redwood ~ 1, ctrl=list(rmin=RMIN,q=1/2))
  if(fit$Fit$mcfit$ctrl$rmin != RMIN)
    stop("kppm did not handle parameter 'rmin' in argument 'ctrl' ")
  fit <- kppm(redwood ~ 1, ctrl=list(rmin=0,q=1/2), rmin=RMIN)
  if(fit$Fit$mcfit$ctrl$rmin != RMIN)
    stop("kppm did not handle parameter 'rmin' in argument 'ctrl'")

  RMIN <- 2
  fit <- dppm(swedishpines~1, dppGauss(), ctrl=list(rmin=RMIN,q=1))
  if(fit$Fit$mcfit$ctrl$rmin != RMIN)
    stop("dppm did not handle parameter 'rmin' in argument 'ctrl'")
  fit <- dppm(swedishpines~1, dppGauss(), ctrl=list(rmin=0,q=1), rmin=RMIN)
  if(fit$Fit$mcfit$ctrl$rmin != RMIN)
    stop("dppm did not handle argument 'rmin'")
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
  ## unsupported options that give a warning
  spatstat.options(kppm.canonical=TRUE, kppm.adjusted=TRUE)
  futXX1 <- kppm(redwood, clusters="MatClust") 
  futXX2 <- kppm(redwood, clusters="MatClust", method="palm")
  futXX3 <- kppm(redwood, clusters="MatClust", method="clik2")
  jpines <- residualspaper$Fig1
  fut <- dppm(jpines ~ 1, dppGauss)
  print(fut)
  spatstat.options(kppm.canonical=FALSE, kppm.adjusted=FALSE)
})

local({
  #' cover a few code blocks
  fut <- kppm(redwood ~ x, method="clik")
  print(summary(fut))
  a <- residuals(fut)
  fut2 <- kppm(redwood ~ x, "LGCP", method="palm")
  print(summary(fut2))
  b <- residuals(fut2)
  #'
  po <- ppm(redwood ~ 1)
  A <- kppmComLik(redwood, Xname="redwood", po=po, clusters="Thomas",
                  statistic="pcf", statargs=list(), control=list(),
                  weightfun=NULL, rmax=0.1)
  A <- kppmPalmLik(redwood, Xname="redwood", po=po, clusters="Thomas",
                   statistic="pcf", statargs=list(), control=list(),
                   weightfun=NULL, rmax=0.1)
})

reset.spatstat.options()
  
