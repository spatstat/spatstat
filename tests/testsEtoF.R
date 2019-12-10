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
#  $Revision: 1.20 $  $Date: 2019/12/02 06:29:20 $
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

#' check handling of 'dangerous' cases
fut <- ppm(redwood ~ x)
Ek <- envelope(fut, Kinhom, update=FALSE, nsim=4)
kfut <- kppm(redwood3 ~ x)
Ekk <- envelope(kfut, Kinhom, lambda=density(redwood3), nsim=7)

# check conditional simulation
e1 <- envelope(cells, Kest, nsim=4, fix.n=TRUE)
e2 <- envelope(amacrine, Kest, nsim=4, fix.n=TRUE)
e3 <- envelope(amacrine, Kcross, nsim=4, fix.marks=TRUE)
e4 <- envelope(finpines, Kest, nsim=4, fix.n=TRUE) # multiple columns of marks
e5 <- envelope(finpines, Kest, nsim=4, fix.marks=TRUE)
fit <- ppm(japanesepines ~ 1, Strauss(0.04))
e6 <- envelope(fit, Kest, nsim=4, fix.n=TRUE)
fit2 <- ppm(amacrine ~ 1, Strauss(0.03))
e7 <- envelope(fit2, Gcross, nsim=4, fix.marks=TRUE)

# check pooling of envelopes in global case
E1 <- envelope(cells, Kest, nsim=5, savefuns=TRUE, global=TRUE)
E2 <- envelope(cells, Kest, nsim=12, savefuns=TRUE, global=TRUE)
p12 <- pool(E1, E2)
p12 <- pool(E1, E2, savefuns=TRUE)
F1 <- envelope(cells, Kest, nsim=5,
               savefuns=TRUE, savepatterns=TRUE, global=TRUE)
F2 <- envelope(cells, Kest, nsim=12,
               savefuns=TRUE, savepatterns=TRUE, global=TRUE)
p12 <- pool(F1, F2)
p12 <- pool(F1, F2, savefuns=TRUE, savepatterns=TRUE)
E1r <- envelope(cells, Kest, nsim=5, savefuns=TRUE, global=TRUE,
                ginterval=c(0.05, 0.15))
E2r <- envelope(cells, Kest, nsim=12, savefuns=TRUE, global=TRUE,
                ginterval=c(0.05, 0.15))
p12r <- pool(E1r, E2r)
})

local({
  #' as.data.frame.envelope
  Nsim <- 5
  E <- envelope(cells, nsim=Nsim, savefuns=TRUE)
  A <- as.data.frame(E)
  B <- as.data.frame(E, simfuns=TRUE)
  stopifnot(ncol(B) - ncol(A) == Nsim)
})

local({
  #' cases not covered elsewhere
  A <- envelope(cells, nsim=5, alternative="less",
                do.pwrong=TRUE, use.theory=FALSE,
                savepatterns=TRUE, savefuns=TRUE)
  print(A)
  B <- envelope(A, nsim=5, savefuns=TRUE)
  D <- envelope(cells, "Lest", nsim=5)
  
  UU <- envelope(cells, nsim=5, foreignclass="ppp", clipdata=TRUE)
  
  AA <- envelope(cells, nsim=5, jsim=5, alternative="greater", global=TRUE)
  AA <- envelope(cells, nsim=5, jsim=5, alternative="less", global=TRUE)
  AA <- envelope(cells, nsim=5, jsim=5, alternative="greater", VARIANCE=TRUE)
  AA <- envelope(cells, nsim=5, jsim=5, alternative="greater", VARIANCE=TRUE)

  fit <- ppm(cells ~ 1, Strauss(0.07))
  U <- envelope(fit, nsim=3, simulate=expression(runifpoint(20)))
  kfit <- kppm(redwood3 ~ x)
  UU <- envelope(kfit, nsim=7, simulate=expression(simulate(kfit, drop=TRUE)))
  VV <- envelope(kfit, nsim=7, weights=1:7)
  MM <- envelope(kfit, nsim=7, Kinhom, lambda=density(redwood3))
  #'  envelopes based on sample variance
  E <- envelope(cells, nsim=8, VARIANCE=TRUE)
  G <- envelope(cells, nsim=8, VARIANCE=TRUE,
                use.theory=FALSE, do.pwrong=TRUE)
  print(G)
  #' summary method
  summary(E)
  summary(envelope(cells, nsim=5, simulate=expression(runifpoint(42))))
  #' weights argument
  H1 <- envelope(cells, nsim=4, weights=npoints, savefuns=TRUE)
  H2 <- envelope(cells, nsim=4, weights=npoints, savefuns=TRUE)
  J1 <- envelope(cells, nsim=4, weights=npoints, VARIANCE=TRUE)
  J2 <- envelope(cells, nsim=4, weights=npoints, VARIANCE=TRUE)
  #' pooling with weights
  H <- pool(H1, H2)
  J <- pool(J1, J2)
  #' pooling envelopes with non-identical attributes
  H0 <- envelope(cells, nsim=4, savefuns=TRUE)
  HH <- pool(H0, H1)
  #' undocumented/secret
  K <- envelope(cells, nsim=4, saveresultof=npoints, collectrubbish=TRUE)
  #' so secret I've even forgotten how to do it
  M <- envelope(cells, nsim=4, internal=list(eject="patterns"))
})

local({
  #' envelope computations in other functions
  P <- lurking(cells, expression(x), envelope=TRUE, nsim=9)
  print(P)
  #' re-using envelope objects in other functions
  A <- envelope(cells, nsim=9, savepatterns=TRUE, savefuns=TRUE)
  S <- lurking(cells, expression(x), envelope=A, nsim=9)
  #' envelope.envelope
  B <- envelope(cells, nsim=5, savepatterns=TRUE, savefuns=FALSE)
  envelope(B)
})

local({
  X <- runiflpp(10, simplenet)
  Xr <- X %mark% runif(10)
  Xc <- X %mark% factor(letters[c(1:4,3,2,4:1)])
  X2 <- X %mark% data.frame(height=runif(10), width=runif(10))

  E  <- envelope(X,  linearK, nsim=9)
  Er <- envelope(Xr, linearK, nsim=9)
  Ec <- envelope(Xc, linearK, nsim=9)
  E2 <- envelope(X2, linearK, nsim=9)
  
  Erf <- envelope(Xr, linearK, nsim=9, fix.n=TRUE)
  E2f <- envelope(X2, linearK, nsim=9, fix.n=TRUE)

  Ecf <- envelope(Xc, linearK,      nsim=9, fix.n=TRUE)
  Ecm <- envelope(Xc, linearKcross, nsim=9, fix.n=TRUE, fix.marks=TRUE)

  fut <- lppm(Xc ~ marks)
  EEf <- envelope(fut, linearK,      fix.n=TRUE)
  EEm <- envelope(fut, linearKcross, fix.n=TRUE, fix.marks=TRUE)
})

local({
  #' Test robustness of envelope() sorting procedure when NA's are present
  #' Fails with spatstat.utils 1.12-0
  set.seed(42)
  EP <- envelope(longleaf, pcf, nsim=10, nrank=2)

  #' Test case when the maximum permitted number of failures is exceeded
  X <- amacrine[1:153] # contains exactly one point with mark='off'
  #' High probability of generating a pattern with no marks = 'off'
  E <- envelope(X, Kcross, nsim=39, maxnerr=2, maxerr.action="warn")
  A <- alltypes(X, Kcross, envelope=TRUE, nsim=39, maxnerr=2)
})

local({
  #' Internals: envelope.matrix
  Y <- matrix(rnorm(200), 10, 20)
  rr <- 1:10
  oo <- rnorm(10)
  zz <- numeric(10)
  E <- envelope(Y, rvals=rr, observed=oo, nsim=10)
  E <- envelope(Y, rvals=rr, observed=oo, jsim=1:10)
  E <- envelope(Y, rvals=rr, observed=oo, theory=zz,
                type="global", use.theory=TRUE)
  E <- envelope(Y, rvals=rr, observed=oo, theory=zz,
                type="global", use.theory=TRUE, nsim=10)
  E <- envelope(Y, rvals=rr, observed=oo, type="global",
                nsim=10, nsim2=10)
  E <- envelope(Y, rvals=rr, observed=oo, type="global",
                jsim=1:10, jsim.mean=11:20)
  print(E)
  E <- envelope(Y, rvals=rr, observed=oo, type="global",
                nsim=10, jsim.mean=11:20)
  E <- envelope(Y, rvals=rr, observed=oo, type="global",
                jsim=1:10, nsim2=10)
})

local({
  #' quirk with handmade summary functions ('conserve' attribute)
  Kdif <- function(X, r=NULL) { # note no ellipsis
    Y <- split(X)
    K1 <- Kest(Y[[1]], r=r)
    K2 <- Kest(Y[[2]], r=r)
    D <- eval.fv(K1-K2)
    return(D)
  }
  envelope(amacrine, Kdif, nsim=3)
})
#'  tests/enveltest.R
#'     Envelope tests (dclf.test, mad.test)
#'     and two-stage tests (bits.test, dg.test, bits.envelope, dg.envelope)
#' 
#'     $Revision: 1.1 $  $Date: 2019/10/09 06:28:21 $ 
#'
require(spatstat)
local({
  #' handling of NA function values (due to empty point patterns)
  set.seed(1234)
  X <- rThomas(5, 0.05, 10) 
  fit <- kppm(X ~ 1, "Thomas")
  set.seed(100000)
  dclf.test(fit)
  set.seed(909)
  dg.test(fit, nsim=9)
})
#
#    tests/factorbugs.R
#
# check for various bugs related to factor conversions
#
#    $Revision: 1.4 $  $Date: 2019/02/10 07:21:02 $
#
require(spatstat)
local({
  ## make a factor image
  m <- factor(rep(letters[1:4], 4))
  Z <- im(m, xcol=1:4, yrow=1:4)
  ## make a point pattern
  set.seed(42)
  X <- runifpoint(20, win=as.owin(Z))
  ## look up the image at the points of X
  ## (a) internal
  ans1 <- lookup.im(Z, X$x, X$y)
  stopifnot(is.factor(ans1))
  ## (b) user level
  ans2 <- Z[X]
  stopifnot(is.factor(ans2))
  ## (c) turn the image into a tessellation
  ##  and apply quadratcount
  V <- tess(image = Z)
  quadratcount(X, tess=V)
  ## (d) pad image
  Y <- padimage(Z, factor("b", levels=levels(Z)))
  stopifnot(Y$type == "factor")
  U <- padimage(Z, "b")
  stopifnot(U$type == "factor")
  ## (e) manipulate levels
  Zb <- relevel(Z, "b")
  Zv <- mergeLevels(Z, vowel="a", consonant=c("b","c","d"))
  P <- X %mark% Z[X]
  Pv <- mergeLevels(P, vowel="a", consonant=c("b","c","d"))
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
#       and options not tested elsewhere
#
#   $Revision: 1.4 $   $Date: 2018/04/13 04:01:25 $
#
require(spatstat)
local({
  ## fast code
  Kb <- Kest(cells, nlarge=0)
  Ku <- Kest(cells, correction="none")
  Kbu <- Kest(cells, correction=c("none", "border"))
  ## slow code, full set of corrections, sqrt transformation, ratios
  Ldd <- Lest(unmark(demopat), correction="all", var.approx=TRUE, ratio=TRUE)
  ## Lotwick-Silverman var approx (rectangular window)
  Loo <- Lest(cells, correction="all", var.approx=TRUE, ratio=TRUE)
  ## Code for large dataset
  nbig <- .Machine$integer.max
  if(!is.null(nbig)) {
    nn <- ceiling(sqrt(nbig))
    if(nn < 1e6) Kbig <- Kest(runifpoint(nn),
                              correction=c("border", "bord.modif", "none"),
                              ratio=TRUE)
  }
  
  ## Kinhom
  lam <- density(cells, at="points", leaveoneout=TRUE)
  ## fast code
  Kib <- Kinhom(cells, lam, nlarge=0)
  Kiu <- Kest(cells, lam, correction="none")
  Kibu <- Kest(cells, lam, correction=c("none", "border"))
  ## slow code
  Lidd <- Linhom(unmark(demopat), sigma=bw.scott)
})


#' tests/formuli.R
#'
#'  Test machinery for manipulating formulae
#' 
#' $Revision: 1.5 $  $Date: 2019/10/05 11:03:26 $

require(spatstat)
local({

  ff <- function(A, deletevar, B) {
    D <- reduceformula(A, deletevar)
    if(!spatstat.utils::identical.formulae(D, B)) {
      AD <- as.expression(substitute(reduceformula(A,d),
                                     list(A=A, d=deletevar)))
      stop(paste(AD, "\n\tyields ", spatstat.utils::pasteFormula(D),
                 " instead of ", spatstat.utils::pasteFormula(B)),
           call.=FALSE)
    }
    invisible(NULL)
  }

  ff(~ x + z, "x", ~z)

  ff(y ~ x + z, "x", y~z)

  ff(~ I(x^2) + z, "x",  ~z)

  ff(y ~ poly(x,2) + poly(z,3), "x", y ~poly(z,3))

  illegal.iformula(~str*g, itags="str", dfvarnames=c("marks", "g", "x", "y"))
})



#
#  tests/func.R
#
#   $Revision: 1.4 $   $Date: 2019/01/14 07:05:27 $
#
#  Tests of 'funxy' infrastructure etc

require(spatstat)
local({
  ## Check the peculiar function-building code in funxy
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
  ## check coordinate extraction from objects
  X <- runifpoint(9)
  Y <- runiflpp(5, simplenet)
  Q <- quadscheme(X)
  a <- F1a(X)
  b <- F1a(Y)
  d <- F1a(Q)
})




##  
##     tests/funnymarks.R
##
## tests involving strange mark values
## $Revision: 1.5 $ $Date: 2018/06/14 01:50:19 $

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
  md <- markformat(endoftime)
  
  ## mark formats
  Z <- Y
  marks(Z) <- marks(Z)[1,,drop=FALSE]
  ms <- markformat(solist(cells, redwood))
  marks(Z) <- factor(1:npoints(Z))
  marks(Z)[12] <- NA
  mz <- is.multitype(Z)
  cZ <- coerce.marks.numeric(Z)
  stopifnot(is.multitype(cells %mark% data.frame(a=factor(1:npoints(cells)))))

  a <- numeric.columns(finpines)
  b1 <- numeric.columns(amacrine)
  b2 <- coerce.marks.numeric(amacrine)
  d <- numeric.columns(cells)
  f <- numeric.columns(longleaf)
})
## 
##    tests/fvproblems.R
##
##    problems with fv and fasp code
##
##    $Revision: 1.10 $  $Date: 2019/12/06 02:51:36 $

require(spatstat)

#' This appears in the workshop notes
#' Problem detected by Martin Bratschi

local({
  Jdif <- function(X, ..., i) {
    Jidot <- Jdot(X, ..., i=i)
    J <- Jest(X, ...)
    dif <- eval.fv(Jidot - J)
    return(dif)
  }
  Z <- Jdif(amacrine, i="on")
})
#'
#'  Test mathlegend code
#'
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

#'
#'  Test quirks related to 'alim' attribute

local({
  K <- Kest(cells)
  attr(K, "alim") <- NULL
  plot(K)
  attr(K, "alim") <- c(0, 0.1)
  plot(tail(K))
})

#'
#' Check that default 'r' vector passes the test for fine spacing

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

##' various functionality in fv.R

local({
  M <- cbind(1:20, matrix(runif(100), 20, 5))
  A <- as.fv(M)
  fvlabels(A) <- c("r","%s(r)", "%s[A](r)", "%s[B](r)", "%s[C](r)", "%s[D](r)")
  A <- rename.fv(A, "M", quote(M(r)))
  A <- tweak.fv.entry(A, "V1", new.tag="r")
  A[,3] <- NULL
  A$hogwash <- runif(nrow(A))
  #' bind.fv with qualitatively different functions
  GK <- harmonise(G=Gest(cells), K=Kest(cells))
  G <- GK$G
  K <- GK$K
  ss <- c(rep(TRUE, nrow(K)-10), rep(FALSE, 10))
  U <- bind.fv(G, K[ss, ], clip=TRUE)
  #'
  H <- rebadge.as.crossfun(K, "H", "inhom", 1, 2)
  H <- rebadge.as.dotfun(K, "H", "inhom", 3)
})

local({
  ## bug in Jmulti.R colliding with breakpts.R
  B <- owin(c(0,3), c(0,10))
  Y <- superimpose(A=runifpoint(1212, B), B=runifpoint(496, B))
  JDX <- Jdot(Y)
  JCX <- Jcross(Y)
  Jdif <- function(X, ..., i) {
    Jidot <- Jdot(X, ..., i=i)
    J <- Jest(X, ...)
    dif <- eval.fv(Jidot - J)
    return(dif)
  }
  E <- envelope(Y, Jdif, nsim=19, i="A", simulate=expression(rlabel(Y)))
})

local({
  #' fasp axes and title
  a <- alltypes(amacrine)
  a$title <- NULL
  plot(a, samex=TRUE, samey=TRUE)
})
