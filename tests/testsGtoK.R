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
#  $Revision: 1.3 $  $Date: 2014/08/25 04:43:07 $
#

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


#
#  tests/imageops.R
#
#   $Revision: 1.7 $   $Date: 2015/12/29 08:54:49 $
#

require(spatstat)
local({
  A <- as.im(owin())
  B <- as.im(owin(c(1.1, 1.9), c(0,1)))
  Z <- imcov(A, B)
  stopifnot(abs(max(Z) - 0.8) < 0.1)

  ## handling images with 1 row or column
  ycov <- function(x, y) y
  E <- as.im(ycov, owin(), dimyx = c(2,1))
  G <- cut(E, 2)
  H <- as.tess(G)

  E12 <- as.im(ycov, owin(), dimyx = c(1,2))
  G12 <- cut(E12, 2)
  H12 <- as.tess(G12)

  ##
  d <- distmap(cells, dimyx=32)
  Z <- connected(d <= 0.06, method="interpreted")
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
})#
# tests/kppm.R
#
# $Revision: 1.15 $ $Date: 2016/09/13 02:30:05 $
#
# Test functionality of kppm that depends on RandomFields
# Test update.kppm for old style kppm objects

require(spatstat)
local({

 fit <- kppm(redwood, ~1, "Thomas")
 fitx <- update(fit, ~ . + x)
 fitM <- update(fit, clusters="MatClust")
 fitC <- update(fit, cells)
 fitCx <- update(fit, cells ~ x)

 # improve.kppm
 fitI <- update(fit, improve.type="quasi")
 fitxI <- update(fitx, improve.type="quasi")
 # vcov.kppm
 vc <- vcov(fitxI)

  # plot.kppm including predict.kppm
 fitMC <- kppm(redwood ~ x, "Thomas")
 fitCL <- kppm(redwood ~ x, "Thomas", method="c")
 fitPA <- kppm(redwood ~ x, "Thomas", method="p")
 plot(fitMC)
 plot(fitCL)
 plot(fitPA)

 # fit with composite likelihood method [thanks to Abdollah Jalilian]
 fut <- kppm(redwood ~ x, "VarGamma", method="clik2", nu.ker=-3/8)
  
 if(require(RandomFields)) {
   fit0 <- kppm(redwood ~1, "LGCP")
   Y0 <- simulate(fit0)[[1]]
   stopifnot(is.ppp(Y0))

   # fit LGCP using K function: slow
   fit1 <- kppm(redwood ~x, "LGCP",
                covmodel=list(model="matern", nu=0.3),
                control=list(maxit=3))
   Y1 <- simulate(fit1)[[1]]
   stopifnot(is.ppp(Y1))

   # fit LGCP using pcf
   fit1p <- kppm(redwood ~x, "LGCP",
                 covmodel=list(model="matern", nu=0.3),
                 statistic="pcf")
   Y1p <- simulate(fit1p)[[1]]
   stopifnot(is.ppp(Y1p))
  
   # ... and Abdollah's code

   fit2 <- kppm(redwood ~x, cluster="Cauchy", statistic="K")
   Y2 <- simulate(fit2)[[1]]
   stopifnot(is.ppp(Y2))

 }

})


