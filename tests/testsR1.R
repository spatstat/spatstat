#'  tests/randoms.R
#'   Further tests of random generation code
#'  $Revision: 1.11 $ $Date: 2019/12/11 00:30:12 $

require(spatstat)
local({
  A <- runifrect(6, nsim=2)
  A <- runifdisc(6, nsim=2)
  A <- runifpoispp(5, nsim=2)
  A <- runifpoispp(0, nsim=2)
  A <- rSSI(0.05, 6, nsim=2)
  A <- rSSI(0.05, 10, win=square(c(-0.5, 1.5)), x.init=A[[1]], nsim=2)  
  A <- rstrat(nx=4, nsim=2)
  A <- rsyst(nx=4, nsim=2)
  A <- rthin(cells, P=0.5, nsim=2)
  A <- rthin(cells, runif(42))
  A <- rthin(cells[FALSE], P=0.5, nsim=2)
  A <- rjitter(cells, nsim=2, retry=FALSE)
  A <- rjitter(cells[FALSE])
  A <- rcell(square(1), nx=5, nsim=2)

  f <- function(x,y) { 10*x }
  Z <- as.im(f, square(1))
  A <- rpoint(n=6, f=f, fmax=10, nsim=2)
  A <- rpoint(n=6, f=Z, fmax=10, nsim=2)
  A <- rpoint(n=0, f=f, fmax=10, nsim=2)
  A <- rpoint(n=0, f=Z, fmax=10, nsim=2)

  op <- spatstat.options(fastpois=FALSE)
  A <- runifpoispp(5, nsim=2)
  A <- rpoispp(Z)
  spatstat.options(op)

  b3 <- box3(c(0,1))
  b4 <- boxx(c(0,1), c(0,1), c(0,1), c(0,1))
  b5 <- c(0, 2, 0, 2)
  X <- rMaternInhibition(2, kappa=20, r=0.1, win=b3)
  Y <- rMaternInhibition(2, kappa=20, r=0.1, win=b4)
  Y <- rMaternInhibition(2, kappa=20, r=0.1, win=b5, nsim=2)

  X <- rSSI(0.05, 6)
  Y <- rSSI(0.05, 6, x.init=X) # no extra points

  Z <- rlabel(finpines)
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
  G2 <- rpoisppOnLines(list(A=f1, B=f2), L, lmax=max(sapply(ZZ, max)))
})

local({
  #' cluster models + bells + whistles
  X <- rThomas(10, 0.2, 5, saveLambda=TRUE)
  if(is.null(attr(X, "Lambda")))
    stop("rThomas did not save Lambda image")
  Y <- rThomas(0, 0.2, 5, saveLambda=TRUE)
  if(is.null(attr(Y, "Lambda")))
    stop("rThomas did not save Lambda image when kappa=0")
  X <- rMatClust(10, 0.05, 4, saveLambda=TRUE)
  X <- rCauchy(30, 0.01, 5, saveLambda=TRUE)
  X <- rVarGamma(30, 2, 0.02, 5, saveLambda=TRUE)
  Z <- as.im(function(x,y){ 5 * exp(2 * x - 1) }, owin())
  Y <- rThomas(10, 0.2, Z, saveLambda=TRUE)
  Y <- rMatClust(10, 0.05, Z, saveLambda=TRUE)
  Y <- rCauchy(30, 0.01, Z, saveLambda=TRUE)
  Y <- rVarGamma(30, 2, 0.02, Z, saveLambda=TRUE)
})

#' perfect simulation code
expandwinPerfect(letterR, 2, 3)

reset.spatstat.options()


#'  tests/resid.R
#'
#'  Stuff related to residuals and residual diagnostics
#'
#'   $Revision: 1.5 $  $Date: 2019/12/21 04:59:50 $
#'

require(spatstat)
local({
  fit <- ppm(cells ~x, Strauss(r=0.15))
  diagnose.ppm(fit, cumulative=FALSE)
  diagnose.ppm(fit, cumulative=FALSE, type="pearson")
  rr <- residuals(fit, quad=quadscheme(cells, nd=128))

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

  #' entangled terms in model
  U <- as.im(1, owin())
  Z <- as.im(function(x,y) x, owin())
  X <- runifpoint(40)
  fut <- ppm(X ~ Z:U)
  a <- parres(fut, "Z")
  futoff <- ppm(X ~ offset(Z*U))
  a <- parres(futoff, "Z")
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
