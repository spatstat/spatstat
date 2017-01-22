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
#  $Revision: 1.8 $  $Date: 2017/01/22 01:50:49 $
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

  # apply different discretisation rules
  Z <- density(cells, 0.05, fractional=TRUE)
  Z <- density(cells, 0.05, preserve=TRUE)
  Z <- density(cells, 0.05, fractional=TRUE, preserve=TRUE)
        
  ## compare density.ppp(at="points") results with different algorithms
  crosscheque <- function(expr) {
    e <- as.expression(substitute(expr))
    ename <- sQuote(deparse(substitute(expr)))
    ## interpreted R
    opa <- spatstat.options(densityC=FALSE, densityTransform=FALSE)
    val.interpreted <- eval(e)
    ## established C algorithm 'denspt'
    spatstat.options(densityC=TRUE, densityTransform=FALSE)
    val.C <- eval(e)
    ## new C algorithm 'Gdenspt' using transformed coordinates
    spatstat.options(densityC=TRUE, densityTransform=TRUE)
    val.Transform <- eval(e)
    spatstat.options(opa)
    if(max(abs(val.interpreted - val.C)) > 0.001)
      stop(paste("Numerical discrepancy between R and C algorithms in",
                 ename))
    if(max(abs(val.C - val.Transform)) > 0.001)
      stop(paste("Numerical discrepancy between C algorithms",
                 "using transformed and untransformed coordinates in",
                 ename))
    invisible(NULL)
  }

  crosscheque(density(redwood, at="points", sigma=0.13, edge=FALSE))

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
  Y <- longleaf %mark% runif(npoints(longleaf), min=41, max=43)

  Z <- Smooth(Y, 5, at="points", leaveoneout=TRUE)
  rZ <- range(Z)
  if(rZ[1] < 40 || rZ[2] > 44)
    stop("Implausible results from Smooth.ppp(at=points, leaveoneout=TRUE)")

  Z <- Smooth(Y, 5, at="points", leaveoneout=FALSE)
  rZ <- range(Z)
  if(rZ[1] < 40 || rZ[2] > 44)
    stop("Implausible results from Smooth.ppp(at=points, leaveoneout=FALSE)")

  ## compare Smooth.ppp results with different algorithms
  crosscheque(Smooth(longleaf, at="points", sigma=6))

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
  # from Cenk Icos
  X <- runifpoint(100, owin(c(0,3), c(0,10)))
  FX <- Fest(X)
  FXr <- Fest(X, r=FX$r)
  JX <- Jest(X)
})

  
