## badwindowcheck.R
## $Revision: 1.2 $  $Date: 2014/01/27 07:18:41 $
##

require(spatstat)
local({
  ## Simple example of self-crossing polygon
  x <- read.table("selfcross.txt", header=TRUE)
  ## Auto-repair
  w <- owin(poly=x)

  ## Real data involving various quirks
  b <- read.table("badwindow.txt", header=TRUE)
  b <- split(b, factor(b$i))
  b <- lapply(b, function(z) { as.list(z[,-3]) })
  ## make owin without checking
  W <- owin(poly=b, check=FALSE)
  ## Apply stringent checks
  owinpolycheck(W,verbose=FALSE)
  ## Auto-repair
  W2 <- owin(poly=b)
})




## tests/cdf.test.R
## check cdf.test with strange data
require(spatstat)
local({
  # Marked point patterns with some marks not represented
  AC <- split(ants, un=FALSE)$Cataglyphis
  AM <- split(ants, un=FALSE)$Messor
  DM <- distmap(AM)
  # should produce a warning, rather than a crash:
  cdf.test(AC, DM)
  # should be OK:
  cdf.test(unmark(AC), DM)
  cdf.test(unmark(AC), DM, "cvm")
  cdf.test(unmark(AC), DM, "ad")
  # linear networks
  set.seed(42)
  X <- runiflpp(20, simplenet)
  fit <- lppm(X ~1)
  cdf.test(fit, "y")
  cdf.test(fit, "y", "cvm")
  cdf.test(fit, "y", "ad")
})

##  tests/closeshave.R
## check 'closepairs/crosspairs' code
## validity and memory allocation
## $Revision: 1.5 $ $Date: 2016/03/28 04:21:07 $

local({
  r <- 0.12
  close.all <- closepairs(redwood, r)
  close.ij <- closepairs(redwood, r, what="indices")
  close.ijd <- closepairs(redwood, r, what="ijd")
  stopifnot(identical(close.ij, close.all[c("i","j")]))
  stopifnot(identical(close.ijd, close.all[c("i","j","d")]))

  Y <- split(amacrine)
  on <- Y$on
  off <- Y$off
  cross.all <- crosspairs(on, off, r)
  cross.ij <- crosspairs(on, off, r, what="indices")
  cross.ijd <- crosspairs(on, off, r, what="ijd")
  stopifnot(identical(cross.ij, cross.all[c("i","j")]))
  stopifnot(identical(cross.ijd, cross.all[c("i","j","d")]))

  # closethresh vs closepairs: EXACT agreement
  thresh <- 0.08
  clt <- closethresh(redwood, r, thresh)
  cl <- with(closepairs(redwood, r),
             list(i=i, j=j, th = (d <= thresh)))
  if(!identical(cl, clt))
    stop("closepairs and closethresh disagree")

  # compare with older, slower code
  reordered <- function(a) {
    o <- with(a, order(i,j))
    as.list(as.data.frame(a)[o,,drop=FALSE])
  }
  samesame <- function(a, b) {
    identical(reordered(a), reordered(b))
  }
  spatstat.options(closepairs.newcode=FALSE)
  old.close.ij <- closepairs(redwood, r, what="indices")
  old.cross.ij <- crosspairs(on, off, r, what="indices")
  stopifnot(samesame(close.ij, old.close.ij))
  stopifnot(samesame(cross.ij, old.cross.ij))
  spatstat.options(closepairs.newcode=TRUE)
  
  # Rasmus' example
  R <- 0.04
  U <- as.ppp(gridcenters(owin(), 50, 50), W=owin())
  cp <- crosspairs(U, U, R)
  G <- matrix(0, npoints(U), npoints(U))
  G[cbind(cp$i, cp$j)] <- 1
  if(!isSymmetric(G))
    stop("crosspairs is not symmetric in Rasmus example")

})
## tests/colour.R
##
## $Revision: 1.1 $ $Date: 2015/12/29 08:54:49 $
##

require(spatstat)

local({
   f <- function(n) grey(seq(0,1,length=n))
   z <- to.grey(f)
})
# tests/correctC.R
# check for agreement between C and interpreted code
# for interpoint distances etc.
# $Revision: 1.4 $ $Date: 2015/12/29 08:54:49 $

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

  # whist
  set.seed(98123)
  x <- runif(1000)
  w <- sample(1:5, 1000, replace=TRUE)
  b <- seq(0,1,length=101)
  op <- spatstat.options(Cwhist=TRUE)
  aT <- whist(x,b,w)
  spatstat.options(Cwhist=FALSE)
  aF <- whist(x,b,w)
  if(!all(aT == aF))
    stop("Algorithms for whist disagree")
  spatstat.options(op)
})
#
#  tests/density.R
#
#  Test behaviour of density methods and inhomogeneous summary functions
#
#  $Revision: 1.6 $  $Date: 2016/07/06 03:42:14 $
#

require(spatstat)

local({

  # test all cases of density.ppp
  
  tryit <- function(...) {
    Z <- density(cells, ..., at="pixels")
    Z <- density(cells, ..., at="points")
    return(invisible(NULL))
  }
  
  tryit(0.05)
  tryit(0.05, diggle=TRUE)
  tryit(0.05, se=TRUE)
  tryit(varcov=diag(c(0.05^2, 0.07^2)))
  tryit(0.05, weights=data.frame(a=1:42,b=42:1))
  tryit(0.05, weights=expression(x))

  ## compare results with different algorithms
  opa <- spatstat.options(densityC=FALSE, densityTransform=FALSE)
  val.interpreted <- density(redwood, at="points", sigma=0.13, edge=FALSE)
  spatstat.options(densityC=TRUE, densityTransform=FALSE)
  val.C <- density(redwood, at="points", sigma=0.13, edge=FALSE)
  spatstat.options(densityC=TRUE, densityTransform=TRUE)
  val.Transform <- density(redwood, at="points", sigma=0.13, edge=FALSE)
  spatstat.options(opa)

  if(max(abs(val.interpreted - val.C)) > 0.001)
    stop(paste("Numerical discrepancy between R and C algorithms",
               "in density.ppp(at=points)"))
  if(max(abs(val.C - val.Transform)) > 0.001)
    stop(paste("Numerical discrepancy between C algorithms",
               "using transformed and untransformed coordinates",
               "in density.ppp(at=points)"))
  
  lam <- density(redwood)
  K <- Kinhom(redwood, lam)
  
  lamX <- density(redwood, at="points")
  KX <- Kinhom(redwood, lamX)

  ## test all code cases of new 'relrisk.ppp' algorithm
  pants <- function(..., X=ants) {
    a <- relrisk(X, sigma=100, se=TRUE, ...)
    return(TRUE)
  }
  pants()
  pants(casecontrol=FALSE)
  pants(relative=TRUE)
  pants(casecontrol=FALSE, relative=TRUE)
  pants(at="points")
  pants(casecontrol=FALSE,at="points")
  pants(relative=TRUE,at="points")
  pants(casecontrol=FALSE, relative=TRUE,at="points")

  ## more than 2 types
  pants(X=sporophores)
  pants(X=sporophores, at="points")
  pants(X=sporophores, relative=TRUE, at="points")

  ## Smooth.ppp
  Z <- Smooth(longleaf, 5, diggle=TRUE)
  Z <- Smooth(longleaf, 1e-6) # generates warning about small bandwidth

  ## Smooth.ppp(at='points')
  X <- longleaf %mark% 42
  Y <- longleaf %mark% runif(npoints(longleaf), min=41, max=43)

  ZX <- Smooth(X, 5, at="points", leaveoneout=TRUE)
  ZY <- Smooth(Y, 5, at="points", leaveoneout=TRUE)
  rZY <- range(ZY)
  if(rZY[1] < 40 || rZY[2] > 44)
    stop("Implausible results from Smooth.ppp(at=points, leaveoneout=TRUE)")

  ZX <- Smooth(X, 5, at="points", leaveoneout=FALSE)
  ZY <- Smooth(Y, 5, at="points", leaveoneout=FALSE)
  rZY <- range(ZY)
  if(rZY[1] < 40 || rZY[2] > 44)
    stop("Implausible results from Smooth.ppp(at=points, leaveoneout=FALSE)")
  
  ## drop-dimension coding errors
  X <- longleaf
  marks(X) <- cbind(marks(X), 1)
  Z <- Smooth(X, 5)

  ZZ <- bw.smoothppp(finpines, hmin=0.01, hmax=0.012, nh=2) # reshaping problem
})
#'
#'  tests/discarea.R
#'
#'   $Revision: 1.1 $ $Date: 2016/03/28 09:16:03 $
#'

require(spatstat)
local({
  u <- c(0.5,0.5)
  B <- owin(poly=list(x=c(0.3, 0.5, 0.7, 0.4), y=c(0.3, 0.3, 0.6, 0.8)))
  areaGain(u, cells, 0.1, exact=TRUE)
  areaGain(u, cells, 0.1, W=NULL)
  areaGain(u, cells, 0.1, W=B)

  areaLoss(cells[square(0.4)], 0.1, exact=TRUE)
})
#  tests/emptymarks.R
#
# test cases where there are no (rows or columns of) marks
#
#  $Revision: 1.3 $ $Date: 2015/12/29 08:54:49 $

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
#
#  tests/envelopes.R
#
#  Test validity of envelope data
#
#  $Revision: 1.5 $  $Date: 2015/12/29 08:54:49 $
#

require(spatstat)

local({
checktheo <- function(fit) {
  fitname <- deparse(substitute(fit))
  en <- envelope(fit, nsim=4, verbose=FALSE, nrep=1e3)
  nama <- names(en)
  expecttheo <- is.poisson(fit) && is.stationary(fit)
  context <- paste("Envelope of", fitname)
  if(expecttheo) {
    if(!("theo" %in% nama))
      stop(paste(context, "did not contain", sQuote("theo")))
    if("mmean" %in% nama)
      stop(paste(context, "unexpectedly contained", sQuote("mmean")))
  } else {
    if("theo" %in% nama)
      stop(paste(context, "unexpectedly contained", sQuote("theo")))
    if(!("mmean" %in% nama))
      stop(paste(context, "did not contain", sQuote("mmean")))
  }
  cat(paste(context, "has correct format\n"))
}
  
checktheo(ppm(cells))
checktheo(ppm(cells ~x))
checktheo(ppm(cells ~1, Strauss(0.1)))

# check envelope calls from 'alltypes'
a <- alltypes(demopat, Kcross, nsim=4, envelope=TRUE)
b <- alltypes(demopat, Kcross, nsim=4, envelope=TRUE, global=TRUE)

# check 'transform' idioms
A <- envelope(cells, Kest, nsim=4, transform=expression(. - .x))
B <- envelope(cells, Kest, nsim=4, transform=expression(sqrt(./pi) - .x))

#' check savefuns/savepatterns with global 
fit <- ppm(cells~x)
Ef <- envelope(fit, Kest, nsim=4, savefuns=TRUE, global=TRUE)
Ep <- envelope(fit, Kest, nsim=4, savepatterns=TRUE, global=TRUE)

# check conditional simulation
e1 <- envelope(cells, Kest, nsim=4, fix.n=TRUE)
e2 <- envelope(amacrine, Kest, nsim=4, fix.n=TRUE)
e3 <- envelope(amacrine, Kcross, nsim=4, fix.marks=TRUE)
fit <- ppm(japanesepines ~ 1, Strauss(0.04))
e4 <- envelope(fit, Kest, nsim=4, fix.n=TRUE)
fit2 <- ppm(amacrine ~ 1, Strauss(0.03))
e5 <- envelope(fit2, Gcross, nsim=4, fix.marks=TRUE)

# check pooling of envelopes in global case
E1 <- envelope(cells, Kest, nsim=5, savefuns=TRUE, global=TRUE)
E2 <- envelope(cells, Kest, nsim=12, savefuns=TRUE, global=TRUE)
p12 <- pool(E1, E2)
E1r <- envelope(cells, Kest, nsim=5, savefuns=TRUE, global=TRUE,
                ginterval=c(0.05, 0.15))
E2r <- envelope(cells, Kest, nsim=12, savefuns=TRUE, global=TRUE,
                ginterval=c(0.05, 0.15))
p12r <- pool(E1r, E2r)
})
#
#    tests/factorbugs.R
#
# check for various bugs related to factor conversions
#
#    $Revision: 1.3 $  $Date: 2015/12/29 08:54:49 $
#
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
  # (d) pad image
  Y <- padimage(Z, factor("b", levels=levels(Z)))
  stopifnot(Y$type == "factor")
  U <- padimage(Z, "b")
  stopifnot(U$type == "factor")
})


#
#    tests/fastgeyer.R
#
# checks validity of fast C implementation of Geyer interaction
#
#    $Revision: 1.3 $  $Date: 2015/12/29 08:54:49 $
#
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
# and again for sat < 1
# (spotted by Rolf)  
  G <- Geyer(0.11, 0.5)
  a <- pairsat.family$eval(X,U,EP,G$pot,G$par,"border")
  b <-          G$fasteval(X,U,EP,G$pot,G$par,"border")
  if(!all(a==b))
    stop("Results of Geyer()$fasteval and pairsat.family$eval do not match when sat < 1")
})

#
#  tests/fastK.R
#
# check fast and slow code for Kest
#
#   $Revision: 1.2 $   $Date: 2015/12/29 08:54:49 $
#
require(spatstat)
local({
  ## fast code
  Kb <- Kest(cells, nlarge=0)
  Ku <- Kest(cells, correction="none")
  Kbu <- Kest(cells, correction=c("none", "border"))
  ## slow code, full set of corrections, sqrt transformation
  Ldd <- Lest(unmark(demopat), correction="all", var.approx=TRUE)
})


#' tests/formuli.R
#'
#'  Test machinery for manipulating formulae
#' 
#' $Revision: 1.2 $  $Date: 2016/03/05 02:24:32 $

require(spatstat)
local({

  ff <- function(A, deletevar, B) {
    D <- reduceformula(A, deletevar)
    if(!identical.formulae(D, B)) {
      AD <- as.expression(substitute(reduceformula(A,d),
                                     list(A=A, d=deletevar)))
      stop(paste(AD, "\n\tyields ", pasteFormula(D),
                 " instead of ", pasteFormula(B)),
           call.=FALSE)
    }
    invisible(NULL)
  }

  ff(~ x + z, "x", ~z)

  ff(y ~ x + z, "x", y~z)

  ff(~ I(x^2) + z, "x",  ~z)

  ff(y ~ poly(x,2) + poly(z,3), "x", y ~poly(z,3))

})



#
#  tests/func.R
#
#   $Revision: 1.3 $   $Date: 2016/06/10 15:04:08 $
#
#  Tests of 'funxy' infrastructure etc

require(spatstat)
local({
  W <- square(1)
  f1a <- function(x, y) sqrt(x^2 + y^2)
  f1b <- function(x, y) { sqrt(x^2 + y^2) }
  f2a <- function(x, y) sin(x)
  f2b <- function(x, y) { sin(x) } 
  f3a <- function(x, y) sin(x) + cos(x) 
  f3b <- function(x, y) { sin(x) + cos(x) } 
  f4a <- function(x, y) { z <- x + y ; z }
  f4b <- function(x, y) { x + y } 
  F1a <- funxy(f1a, W)
  F1b <- funxy(f1b, W)
  F2a <- funxy(f2a, W)
  F2b <- funxy(f2b, W)
  F3a <- funxy(f3a, W)
  F3b <- funxy(f3b, W)
  F4a <- funxy(f4a, W)
  F4b <- funxy(f4b, W)
  stopifnot(identical(F1a(cells), F1b(cells)))
  stopifnot(identical(F2a(cells), F2b(cells)))
  stopifnot(identical(F3a(cells), F3b(cells)))
  stopifnot(identical(F4a(cells), F4b(cells)))
})




##  
##     tests/funnymarks.R
##
## tests involving strange mark values
## $Revision: 1.3 $ $Date: 2015/12/29 08:54:49 $

require(spatstat)
local({
  ## ppm() where mark levels contain illegal characters
  hyphenated <- c("a", "not-a")
  spaced <- c("U", "non U")
  suffixed <- c("a+", "a*")
  charred <- c("+", "*")

  irad <- matrix(0.1, 2,2)
  hrad <- matrix(0.005, 2, 2)

  tryit <- function(types, X, irad, hrad) { 
    levels(marks(X)) <- types
    fit <- ppm(X ~marks + polynom(x,y,2),
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

  ## marks which are dates
  X <- cells
  n <- npoints(X)
  endoftime <- rep(ISOdate(2001,1,1), n)
  eotDate   <- rep(as.Date("2001-01-01"), n)
  markformat(endoftime)
  markformat(eotDate)
  marks(X) <- endoftime
  print(X)
  Y <- X %mark% data.frame(id=1:42, date=endoftime, dd=eotDate)
  print(Y)
})
# 
#    tests/fvproblems.R
#
#    $Revision: 1.7 $  $Date: 2016/03/08 00:26:23 $
#

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

  levels(marks(amacrine)) <- c("Nasticreechia krorluppia", "Homo habilis")
  plot(Kcross(amacrine))
  plot(alltypes(amacrine, "K"))
  plot(alltypes(amacrine, "J"))
  plot(alltypes(amacrine, pcfcross))
})

#
#  Test quirks related to 'alim' attribute

local({
  K <- Kest(cells)
  attr(K, "alim") <- NULL
  plot(K)
  attr(K, "alim") <- c(0, 0.1)
  plot(tail(K))
})

#
# Check that default 'r' vector passes the test for fine spacing

local({
  a <- Fest(cells)
  A <- Fest(cells, r=a$r)
  b <- Hest(heather$coarse)
  B <- Hest(heather$coarse, r=b$r)
})

  
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




#
# tests/kppm.R
#
# $Revision: 1.11 $ $Date: 2016/03/04 10:48:03 $
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

 if(require(RandomFields) && RandomFieldsSafe()) {

    fit0 <- kppm(redwood ~1, "LGCP")
    Y0 <- simulate(fit0)[[1]]
    stopifnot(is.ppp(Y0))
    
    fit1 <- kppm(redwood ~x, "LGCP",
                covmodel=list(model="matern", nu=0.3),
                control=list(maxit=5))
    Y1 <- simulate(fit1)[[1]]
    stopifnot(is.ppp(Y1))

# ... and Abdollah's code

    fit2 <- kppm(redwood ~x, cluster="Cauchy", statistic="K")
    Y2 <- simulate(fit2)[[1]]
    stopifnot(is.ppp(Y2))
  }

 # improve.kppm
 fitI <- update(fit, improve.type="quasi")
 fitxI <- update(fitx, improve.type="quasi")
 # vcov.kppm
 vc <- vcov(fitxI)
 
})


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
#  $Revision: 1.6 $  $Date: 2016/02/01 09:44:54 $


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


