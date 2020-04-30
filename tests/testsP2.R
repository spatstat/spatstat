#'
#'   Header for all (concatenated) test files
#'
#'   Require spatstat.
#'   Obtain environment variable controlling tests.
#'
#'   $Revision: 1.5 $ $Date: 2020/04/30 05:31:37 $

require(spatstat)
FULLTEST <- (nchar(Sys.getenv("SPATSTAT_TEST", unset="")) > 0)
ALWAYS   <- TRUE
cat(paste("--------- Executing",
          if(FULLTEST) "** ALL **" else "**RESTRICTED** subset of",
          "test code -----------\n"))
# 
#   tests/ppmBadData.R
#
# $Revision: 1.6 $ $Date: 2020/04/30 05:23:52 $

# Testing robustness of ppm and support functions
# when data are rubbish

local({
  if(ALWAYS) {
    ## from Rolf: very large proportion of data is NA
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
  }

  if(ALWAYS) {
    ## from Andrew Bevan: numerical overflow, ill-conditioned Fisher information
  SEED <- 42
    nongranite<- owin(poly = list(x = c(0, 8500, 7000, 6400, 6400, 6700, 7000, 7200, 7300, 8000, 8100, 8800, 9500, 10000, 10000, 0), y = c(0, 0, 2000, 3800, 4000, 5000, 6500, 7400, 7500, 8000, 8100, 9000, 9500, 9600, 10000, 10000)))
    ## Trend on raster grid
    rain <- as.im(X=function(x,y) { x^2 + y^2 }, W=nongranite, dimyx=100)
    ## Generate a point pattern via a Lennard-Jones process
    set.seed(SEED)
    mod4<- rmhmodel(cif="lennard",
                    par=list(beta=1, sigma=250, epsilon=2.2),
                    trend=rain, w=nongranite)
    ljtr<- rmh(mod4, start=list(n.start=80), control=list(p=1, nrep=1e5))

    ## Fit a point process model to the pattern with rain as a covariate
    ## NOTE INCORRECT TREND FORMULA
    ljtrmod <- ppm(ljtr, trend= ~ Z, interaction=NULL, data=list(Z=rain))
    ss <- summary(ljtrmod)
  }
  
  if(FULLTEST) {
    ## From Ege
    ## Degenerate but non-null argument 'covariates'
    xx <- list()
    names(xx) <- character(0)
    fit <- ppm(cells ~x, covariates = xx)
    st <- summary(fit) 
  }

})



#'   tests/ppmclass.R
#'
#'   Class support for ppm
#'
#'   $Revision: 1.7 $ $Date: 2020/04/30 05:23:52 $

if(FULLTEST) {
local({
  #' (1) print.ppm, summary.ppm, print.summary.ppm
  Z <- as.im(function(x,y){x}, Window(cells))
  fitZ <- ppm(cells ~ Z)
  print(fitZ)
  print(summary(fitZ))
  #' logistic
  fitl <- ppm(swedishpines ~ x+y, method="logi")
  print(fitl)
  print(summary(fitl))
  #' Model with covariate arguments
  f <- function(x,y,b) { x+b }
  fitf <- ppm(cells ~ f, covfunargs=list(b=1))
  print(fitf)
  print(summary(fitf))
  #' Invalid model
  fitN <- ppm(redwood ~ 1, Strauss(0.1))
  print(fitN)
  print(summary(fitN))
  #' standard errors in output
  fat <- ppm(cells ~ x, Strauss(0.12))
  op <- spatstat.options(print.ppm.SE='always')
  print(fat)
  spatstat.options(print.ppm.SE='never')
  print(fat)
  print(fitZ)
  spatstat.options(op)

  ## (2) plot.ppm
  plot(fitZ)
  plot(fat, trend=FALSE, cif=FALSE, se=FALSE)

  ## (3) emend.ppm
  fitZe <- emend(fitZ, trace=TRUE)
  ZZ <- Z
  fitZZ <- ppm(cells ~ Z + ZZ)
  fitZZe <- emend(fitZZ, trace=TRUE)
  fitOK  <- ppm(redwood ~1, Strauss(0.1), emend=TRUE)
  print(fitOK)
  fitNot <- ppm(redwood ~1, Strauss(0.1))
  fitSlow <- emend(fitNot, trace=TRUE)
  print(fitSlow)
  op <- spatstat.options(project.fast=TRUE)
  fitFast <- emend(fitNot, trace=TRUE)
  print(fitFast)
  fitZZe <- emend(fitZZ, trace=TRUE)
  spatstat.options(op)
  
  #' (4) methods for other generics
  logLik(fitZ, absolute=TRUE)
  unitname(fitZ)
  unitname(fat) <- c("metre", "metres")
  is.expandable(fitf)
  fit0 <- update(fitZ, . ~ 1)
  anova(fit0, fitZ, override=TRUE)

  ## example from Robert Aue - handling offsets
  X <- demohyper$Points[[1]]
  GH <- Hybrid(G=Geyer(r=0.1, sat=3), H=Hardcore(0.01))
  fit <- ppm(X ~ 1, GH)
  valid.ppm(fit)
})

reset.spatstat.options()
}
#
#   tests/ppmgam.R
#
#   Test ppm with use.gam=TRUE
#
#   $Revision: 1.4 $  $Date: 2020/04/30 05:23:52 $
#

if(FULLTEST) {
local({
  fit <- ppm(nztrees ~s(x,y), use.gam=TRUE)
  mm <- model.matrix(fit)
  mf <- model.frame(fit)
  v <- vcov(fit)
  prd <- predict(fit)
})
}
#'
#'  tests/ppmlogi.R
#'
#' Tests of ppm(method='logi')
#'    and related code (predict, leverage etc)
#'
#' $Revision: 1.15 $  $Date: 2020/04/30 05:23:52 $
#'

local({
  if(FULLTEST) {
    fit <- ppm(cells ~x, method="logi")
    f <- fitted(fit)
    p <- predict(fit)
    u <- summary(fit)
    fitS <- ppm(cells ~x, Strauss(0.12), method="logi")
    fS <- fitted(fitS)
    pS <- predict(fitS)
    uS <- summary(fitS)
    print(uS)

    plot(leverage(fit))
    plot(influence(fit))
    plot(dfbetas(fit))
    plot(leverage(fitS))
    plot(influence(fitS))
    plot(dfbetas(fitS))
  }

  if(FULLTEST) {
    #' same with hard core - A1 is singular
    fitH <- ppm(cells ~x, Strauss(0.08), method="logi")
    print(fitH)
    fH <- fitted(fitH)
    pH <- predict(fitH)
    uH <- summary(fitH)
    print(uH)
    plot(leverage(fitH))
    plot(influence(fitH))
    plot(dfbetas(fitH))
  }
  
  if(FULLTEST) {
    #' logistic fit to data frame of covariates
    z <- c(rep(TRUE, 5), rep(FALSE, 5))
    df <- data.frame(A=z + 2* runif(10),
                     B=runif(10))
    Y <- quadscheme.logi(runifpoint(5), runifpoint(5))
    fut <- ppm(Y ~ A+B, data=df, method="logi")
    sf <- summary(fut)
    print(sf)
  }

  if(FULLTEST) {
    #' vblogit code, just to check that it runs.
    fee <- ppm(cells ~ x, method="VBlogi", nd=21)
    print(fee)
    summary(fee)
    logLik(fee)
    AIC(fee)
    extractAIC(fee)
    Z <- predict(fee)
    summary(Z)
    print(fee$internal$glmfit) # print.vblogit
  }
})

#
#   tests/ppmmarkorder.R
#
# $Revision: 1.4 $  $Date: 2020/04/30 05:23:52 $
#
# Test that predict.ppm, plot.ppm and plot.fitin
# tolerate marks with levels that are not in alpha order
#
if(ALWAYS) { # locale-dependent?
local({
  X <- amacrine
  levels(marks(X)) <- c("ZZZ", "AAA")
  fit <- ppm(X ~marks, MultiStrauss(c("ZZZ","AAA"), matrix(0.06, 2, 2)))
  aa <- predict(fit, type="trend")
  bb <- predict(fit, type="cif")
  plot(fit)
  plot(fitin(fit))
})
}


#
#   tests/ppmscope.R
#
#   Test things that might corrupt the internal format of ppm objects
#
#   $Revision: 1.6 $  $Date: 2020/04/30 05:23:52 $
#

if(ALWAYS) {  # dependent on R version?
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
}
grep#
#   tests/ppmtricks.R
#
#   Test backdoor exits, hidden options, internals and tricks in ppm
#
#   $Revision: 1.19 $  $Date: 2020/04/30 05:23:52 $
#
local({

  ## (1) skip.border
  if(ALWAYS) { # needed below
    fit <- ppm(cells, ~1, Strauss(0.1), skip.border=TRUE)
  }
  
  ## (2) subset arguments of different kinds
  if(FULLTEST) {
    fut <- ppm(cells ~ x, subset=(x > 0.5))
    fot <- ppm(cells ~ x, subset=(x > 0.5), method="logi")
    W <- owin(c(0.4, 0.8), c(0.2, 0.7))
    fut <- ppm(cells ~ x, subset=W)
    fot <- ppm(cells ~ x, subset=W, method="logi")
    V <- as.im(inside.owin, Window(cells), w=W)
    fet <- ppm(cells ~ x, subset=V)
    fet <- ppm(cells ~ x, subset=V, method="logi")
  }

  ## (3) profilepl -> ppm
  ##     uses 'skip.border' and 'precomputed'
  ##     also tests scoping for covariates
  if(FULLTEST) {
    splants <- split(ants)
    mess    <- splants[["Messor"]]
    cats    <- splants[["Cataglyphis"]]
    ss      <- data.frame(r=seq(60,120,by=20),hc=29/6)
    dM      <- distmap(mess,dimyx=256)
    mungf    <- profilepl(ss, StraussHard, cats ~ dM)
    mungp   <- profilepl(ss, StraussHard, trend=~dM, Q=cats)
  }
  
  ## (4) splitting large quadschemes into blocks
  if(FULLTEST) {
    mop <- spatstat.options(maxmatrix=5000)
    qr <- quadBlockSizes(quadscheme(cells))
    pr <- predict(ppm(cells ~ x, AreaInter(0.05)))
    spatstat.options(mop)
    qr <- quadBlockSizes(quadscheme(cells))
  }
  
  ## (5) shortcuts in summary.ppm
  ## and corresponding behaviour of print.summary.ppm
  if(FULLTEST) {
    print(summary(fit, quick=TRUE))
    print(summary(fit, quick="entries"))
    print(summary(fit, quick="no prediction"))
    print(summary(fit, quick="no variances"))
  }

  ## (6) suffstat.R
  if(ALWAYS) {
    fitP <- update(fit, Poisson())
    suffstat.poisson(fitP, cells)
    fit0 <- killinteraction(fit)
    suffstat.poisson(fit0, cells)
  }

  ## (7) various support for class ppm
  if(FULLTEST) {
    fut <- kppm(redwood ~ x)
    A <- quad.ppm(fut)
    Z <- as.im(function(x,y){x}, Window(cells))
    fitZ <- ppm(cells ~ Z)
    U <- getppmOriginalCovariates(fitZ)
  }
  
  ## (8) support for class profilepl
  if(FULLTEST) {
    rr <- data.frame(r=seq(0.05, 0.15, by=0.02))
    ps <- profilepl(rr, Strauss, cells)
    ##  plot(ps) ## covered in plot.profilepl.Rd
    simulate(ps, nrep=1e4)
    parameters(ps)
    fitin(ps)
    predict(ps, type="cif")
  }
                    
  ## (9) class 'plotppm'
  if(FULLTEST) {
    fut <- ppm(amacrine ~ marks + polynom(x,y,2), Strauss(0.07))
    p <- plot(fut, plot.it=FALSE)
    print(p)
    plot(p, how="contour")
    plot(p, how="persp")
  }

  ## (10) ppm -> mpl.engine -> mpl.prepare
  if(ALWAYS) { # includes C code
    fit <- ppm(cells, NULL)
    fit <- ppm(cells ~ x, clipwin=square(0.7))
    fit <- ppm(cells ~ x, subset=square(0.7))
    DG <- as.im(function(x,y){x+y < 1}, square(1))
    fit <- ppm(cells ~ x, subset=DG)
    fit <- ppm(cells ~ x, GLM=glm)
    fit <- ppm(cells ~ x, famille=quasi(link='log', variance='mu'))
    fit <- ppm(cells ~ x, Hardcore(0.07), skip.border=TRUE, splitInf=TRUE)
  }
  
  ## (11) unidentifiable model (triggers an error in ppm)
  if(FULLTEST) {
    Q <- quadscheme(cells)
    M <- mpl.prepare(Q, cells, as.ppp(Q), trend=~1, covariates=NULL,
                     interaction=Hardcore(0.3), correction="none")
  }
})

reset.spatstat.options()
#'
#'   tests/ppp.R
#'
#'   $Revision: 1.10 $ $Date: 2020/04/30 05:41:59 $
#'
#'  Untested cases in ppp() or associated code

local({
  X <- runifpoint(10, letterR)
  Y <- runifpoint(10, complement.owin(letterR))

  if(FULLTEST) {
    #' test handling of points out-of-bounds
    df <- rbind(as.data.frame(X), as.data.frame(Y))
    A <- ppp(df$x, df$y, window=letterR, marks=1:20)
    #' test handling of points with bad coordinates
    df$x[1:3] <- c(Inf, NA, NaN)
    df$y[18:20] <- c(Inf, NA, NaN)
    B <- ppp(df$x, df$y, window=letterR, marks=1:20)
    D <- ppp(df$x, df$y, window=letterR, marks=data.frame(id=1:20, u=runif(20)))
  
    #' test print/summary/plot methods on these bad objects
    print(A)
    print(B)
    print(D)
    print(summary(A))
    print(summary(B))
    print(summary(D))
    plot(A)
    plot(B)
    plot(D)
    plot(attr(A, "rejects"))
    plot(attr(B, "rejects"))
    plot(attr(D, "rejects"))
  
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

    #' cut.ppp
    CU <- cut(U, "height")
    CU <- cut(U, breaks=3)

    #' cases of [<-.ppp
    set.seed(999)
    X <- cells
    B <- square(0.2)
    X[B] <- runifpoint(3, B)
    #' checking 'value'
    Y <- flipxy(X)
    X[B] <- Y[square(0.3)]
    ## deprecated use of second argument
    X[,1:4] <- runifpoint(3)  # deprecated
    X[,B] <- runifpoint(3, B) # deprecated 
    X[1:3, B] <- runifpoint(3, B) # deprecated but does not crash
  
    #' test as.ppp for spatial package if it is not installed
    FR <- Frame(letterR)
    as.ppp(list(x=X$x, y=X$y,
                xl=FR$xrange[1], xu=FR$xrange[2],
                yl=FR$yrange[1], yu=FR$yrange[2]))

    #' various utilities
    periodify(cells, 2)
    periodify(demopat, 2)
  }

  if(ALWAYS) {  # involves C code
    a <- multiplicity(finpines)
    a <- multiplicity(longleaf)
  }

  if(FULLTEST) {
    ## superimpose.ppp, extra cases
    X <- runifpoint(20)
    A <- superimpose(cells, X, W="convex")
    A <- superimpose(cells, X, W=ripras)
    B <- superimpose(concatxy(cells), concatxy(X), W=NULL)
    ## superimpose.splitppp
    Y <- superimpose(split(amacrine))

    ## catch outdated usage of scanpp
    d <- system.file("rawdata", "amacrine", package="spatstat.data")
    if(nzchar(d)) {
      W <- owin(c(0, 1060/662), c(0, 1))
      Y <- scanpp("amacrine.txt", dir=d, window=W, multitype=TRUE)
      print(Y)
    }
    ## (bad) usage of cobble.xy
    xx <- runif(10)
    yy <- runif(10)
    W1 <- cobble.xy(xx, yy)
    W2 <- cobble.xy(xx, yy, boundingbox)
    Wnope <- cobble.xy(xx, yy, function(x,y) {cbind(x,y)}, fatal=FALSE)
  }
})

#
# tests/ppx.R
#
# Test operations for ppx objects
#
#  $Revision: 1.6 $ $Date: 2020/04/30 05:41:59 $
#

local({
  if(ALWAYS) {
    ## make data
    df <- data.frame(x=c(1,2,2,1)/4, y=c(1,2,3,1)/4, z=c(2,3,4,3)/5)
    X <- ppx(data=df, coord.type=rep("s", 3), domain=box3())
  }
  if(ALWAYS) {
    #' methods involving C code
    unique(X)
    duplicated(X)
    anyDuplicated(X)
    multiplicity(X)
    uniquemap(X)
  }
  if(FULLTEST) {
    #' general tests
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
  }

  if(FULLTEST) {
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
  }
})
#
# tests/prediction.R
#
# Things that might go wrong with predict()
#
#  $Revision: 1.20 $ $Date: 2020/04/30 05:41:59 $
#

local({
  if(ALWAYS) {
    ## test of 'covfunargs' - platform dependent?
    f <- function(x,y,a){ y - a }
    fit <- ppm(cells ~x + f, covariates=list(f=f), covfunargs=list(a=1/2))
    p <- predict(fit)

    ## prediction involving 0 * NA
    qc <- quadscheme(cells, nd=10)
    r <- minnndist(as.ppp(qc))/10
    fit <- ppm(qc ~ 1, Strauss(r)) # model has NA for interaction coefficient
    p1 <- predict(fit)
    p2 <- predict(fit, type="cif", ngrid=10)
    stopifnot(all(is.finite(as.matrix(p1))))
    stopifnot(all(is.finite(as.matrix(p2))))
  }

  if(FULLTEST) {
    ## test of 'new.coef' mechanism
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
  }

  if(FULLTEST) {
    ## tests of relrisk.ppm
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
    b <- predict(fit, se=TRUE, locations=cells)
    b <- predict(fit, se=TRUE, interval="confidence")
    b <- predict(fit, type="count",                            se=TRUE)
    b <- predict(fit, type="count", window=square(0.5),        se=TRUE)
    b <- predict(fit, type="count", window=quadrats(cells, 3), se=TRUE)
    d <- predict(fit, type="count", interval="prediction", se=TRUE)
    d <- predict(fit, type="count", interval="confidence", se=TRUE)
    d <- predict(fit,               interval="confidence", se=TRUE)
    foot <- ppm(cells ~ x, StraussHard(0.12))
    d <- predict(foot, ignore.hardcore=TRUE)
    dX <- predict(foot, ignore.hardcore=TRUE, locations=cells)
  
    ## superseded usages
    b <- predict(fit, type="se", getoutofjail=TRUE)
    b <- predict(fit, type="se", locations=cells) # warning
    b <- predict(fit, total=TRUE)
    b <- predict(fit, total=square(0.5))
    b <- predict(fit, total=quadrats(cells, 3))

    ## supporting code
    u <- model.se.image(fit, square(0.5))
    u <- model.se.image(fit, square(0.5), what="cv")
    u <- model.se.image(fit, square(0.5), what="ce")
    co <- c(Intercept=5, slope=3, kink=2)
    re <- c("Intercept", "slope")
    a <- fill.coefs(co, re) # warning
    b <- fill.coefs(co, rev(names(co)))
    d <- fill.coefs(co, letters[1:3])
    ## model matrix etc
    v <- model.frame(ppm(cells))
    fut <- ppm(cells ~ x, Strauss(0.1))
    v <- model.matrix(fut, subset=(x<0.5), keepNA=FALSE)
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
    #' given covariate
    dlin <- distfun(copper$SouthLines)
    copfit <- ppm(copper$SouthPoints ~ dlin, Geyer(1,1))
    effectfun(copfit, "dlin")
    effectfun(copfit)
    #' covariate that is not used in model
    effectfun(fut, "y", x=0)
    futS <- ppm(cells ~ 1, Strauss(0.1))
    effectfun(futS, "x")
    effectfun(futS, "y")
    #' factor covariate
    fot <- ppm(amacrine~x+marks)
    effectfun(fot, "marks", x=0.5, se.fit=TRUE)
    #' covariate retained but not used
    W <- Window(swedishpines)
    a <- solist(A=funxy(function(x,y){x < 20}, W),
                B=funxy(function(x,y){factor(x < 20)}, W))
    fvt <- ppm(swedishpines ~ A, data=a, allcovar=TRUE)
    effectfun(fvt, "A",         se.fit=TRUE)
    effectfun(fvt, "B", A=TRUE, se.fit=TRUE)
              
    ## ppm with covariate values in data frame
    X <- rpoispp(42)
    Q <- quadscheme(X)
    weirdfunction <- function(x,y){ 10 * x^2 + 5 * sin(10 * y) }
    Zvalues <- weirdfunction(x.quad(Q), y.quad(Q))
    fot <- ppm(Q ~ y + Z, data=data.frame(Z=Zvalues))
    effectfun(fot, "y", Z=0)
    effectfun(fot, "Z", y=0)

    #' multitype
    modX <- ppm(amacrine ~ polynom(x,2))
    effectfun(modX)
    effectfun(modX, "x")
    modXM <- ppm(amacrine ~ marks*polynom(x,2))
    effectfun(modXM, "x", marks="on")
    modXYM <- ppm(amacrine ~ marks*polynom(x,y,2))
    effectfun(modXYM, "x", y=0, marks="on")
  
    df <- as.data.frame(simulate(modXM, drop=TRUE))
    df$marks <- as.character(df$marks)
    dfpr <- predict(modXM, locations=df)
  }
})

#
#     tests/project.ppm.R
#
#      $Revision: 1.7 $  $Date: 2020/04/30 05:41:59 $
#
#     Tests of projection mechanism
#

local({
  chk <- function(m) {
    if(!valid.ppm(m)) stop("Projected model was still not valid")
    return(invisible(NULL))
  }
  
  if(FULLTEST) {
    ## a very unidentifiable model
    fit <- ppm(cells ~Z, Strauss(1e-06), covariates=list(Z=0))
    chk(emend(fit))
    ## multitype
    r <- matrix(1e-06, 2, 2)
    fit2 <- ppm(amacrine ~1, MultiStrauss(types=c("off", "on"), radii=r))
    chk(emend(fit2))
    ## complicated multitype 
    fit3 <- ppm(amacrine ~1, MultiStraussHard(types=c("off", "on"),
                                              iradii=r, hradii=r/5))
    chk(emend(fit3))

    #' code coverage
    op <- spatstat.options(project.fast=TRUE)
    fut <- emend(fit, trace=TRUE)
    chk(fut)
    spatstat.options(op)

    #' hierarchical
    ra <- r
    r[2,1] <- NA
    fit4 <- ppm(amacrine ~1, HierStrauss(types=c("off", "on"), radii=r))
    chk(emend(fit4))
    #' complicated hierarchical
    fit5 <- ppm(amacrine ~1, HierStraussHard(types=c("off", "on"),
                                             iradii=r, hradii=r/5))
    chk(emend(fit5))
  
    ## hybrids
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
  }
})

reset.spatstat.options()
