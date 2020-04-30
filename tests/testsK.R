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
#'
#'   tests/Kfuns.R
#'
#'   Various K and L functions and pcf
#'
#'   $Revision: 1.39 $  $Date: 2020/04/28 08:18:41 $
#'
#'   Assumes 'EveryStart.R' was run

myfun <- function(x,y){(x+1) * y } # must be outside

local({
  if(FULLTEST) {
    #' supporting code
    rmax.rule("Kscaled", owin(), 42)
    implemented.for.K(c("border", "bord.modif", "translate", "good", "best"),
                      "polygonal", TRUE)
    implemented.for.K(c("border", "bord.modif", "translate", "good", "best"),
                      "mask", TRUE)
    implemented.for.K(c("border", "isotropic"), "mask", TRUE)
    implemented.for.K(c("border", "isotropic"), "mask", FALSE)
    #' shortcuts
    D <- density(cells)
    K <- Kborder.engine(cells, rmax=0.4, weights=D, ratio=TRUE)
    K <- Knone.engine(cells, rmax=0.4, weights=D, ratio=TRUE)
    allcor <- c("none", "border", "bord.modif","isotropic", "translate")
    K <- Krect.engine(cells, rmax=0.4, ratio=TRUE, correction=allcor)
    K <- Krect.engine(cells, rmax=0.4, ratio=TRUE, correction=allcor,
                      weights=D)
    K <- Krect.engine(cells, rmax=0.4, ratio=TRUE, correction=allcor,
                      use.integers=FALSE)
    #' Kest special code blocks
    K <- Kest(cells, var.approx=TRUE, ratio=FALSE)
    Z <- distmap(cells) + 1
    Kb <- Kest(cells, correction=c("border","bord.modif"),
               weights=Z, ratio=TRUE)
    Kn <- Kest(cells, correction="none",
               weights=Z, ratio=TRUE)
    Knb <- Kest(cells, correction=c("border","bord.modif","none"),
                weights=Z, ratio=TRUE)
  }
  if(ALWAYS) {
    bigint <- 50000 # This is only "big" on a 32-bit system where
                    # sqrt(.Machine$integer.max) = 46340.9
    X <- runifpoint(bigint)
    Z <- as.im(1/bigint, owin())
    Kb <- Kest(X, correction=c("border","bord.modif"),
               rmax=0.02, weights=Z, ratio=TRUE)
  }
  if(FULLTEST) {
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
    pv <- pcf(redwood, kernel="rectangular")
    p1 <- pcf(redwood[1])
    #' pcf.fv
    K <- Kest(redwood)
    g <- pcf(K, method="a")
    g <- pcf(K, method="c")
    g <- pcf(K, method="d")
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
  }
  if(ALWAYS) {
    #' edge corrections
    rr <- rep(0.1, npoints(cells))
    eC <- edge.Ripley(cells, rr)
    eI <- edge.Ripley(cells, rr, method="interpreted")
    if(max(abs(eC-eI)) > 0.1)
      stop("Ripley edge correction results do not match")
  }
  if(FULLTEST) {
    a <- rmax.Ripley(square(1))
    a <- rmax.Rigid(square(1))
    a <- rmax.Ripley(as.polygonal(square(1)))
    a <- rmax.Rigid(as.polygonal(square(1)))
    a <- rmax.Ripley(letterR)
    a <- rmax.Rigid(letterR)
  }
  if(ALWAYS) {
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
  }
})

local({
  if(FULLTEST) {
    #' ----  multitype ------
    K <- Kcross(amacrine, correction=c("none", "bord.modif"))
    #' inhomogeneous multitype
    fit <- ppm(amacrine ~ marks)
    K1 <- Kcross.inhom(amacrine, lambdaX=fit)
    K2 <- Kcross.inhom(amacrine, lambdaX=densityfun(amacrine))
    K3 <- Kcross.inhom(amacrine, lambdaX=density(amacrine, at="points"))
    On <- split(amacrine)$on
    Off <- split(amacrine)$off
    K4 <- Kcross.inhom(amacrine, lambdaI=ppm(On), lambdaJ=ppm(Off))
    K5 <- Kcross.inhom(amacrine, correction="bord.modif")
    #' markconnect, markcorr
    M <- markconnect(amacrine, "on", "off", normalise=TRUE)
    M <- markcorr(longleaf, normalise=TRUE,
                  correction=c("isotropic", "translate", "border", "none"))
    M <- markcorr(longleaf, normalise=TRUE, fargs=list())
    #' Kmark (=markcorrint)
    X <- runifpoint(100) %mark% runif(100)
    km <- Kmark(X, f=atan2)
    km <- Kmark(X, f1=sin)
    km <- Kmark(X, f="myfun")
    aa <- Kmark(X, normalise=FALSE, returnL=FALSE)
    aa <- Kmark(X, normalise=FALSE, returnL=TRUE)
    aa <- Kmark(X, normalise=TRUE,  returnL=FALSE)
    aa <- Kmark(X, normalise=TRUE,  returnL=TRUE)
  }
})

local({
  if(FULLTEST) {
    #'    various modified K functions
    #'
    #'   directional K functions
    #'
    a <- Ksector(swedishpines,
                 -pi/2, pi/2, units="radians",
                 correction=c("none", "border", "bord.modif",
                              "Ripley", "translate"),
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
    a <- localLcross(amacrine, from="off", to="off")
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
    d <- densityfun(unmark(amacrine), sigma=0.1)
    dm <- lapply(split(amacrine), densityfun, sigma=0.1)
    h <- resolve.lambda.cross(amacrine, moff, !moff, lambdaX=d)
    h <- resolve.lambda.cross(amacrine, moff, !moff,
                              lambdaI=dm[["off"]], lambdaJ=dm[["on"]])
    h <- resolve.lambda.cross(amacrine, moff, !moff,
                              lambdaX=function(x,y,m){ d(x,y) })
    #'
    #' multitype inhomogeneous pcf
    #'
    g <- pcfcross.inhom(amacrine, 
                        lambdaI=dm[["off"]], lambdaJ=dm[["on"]])

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
    fit <- ppm(cells ~ x, Strauss(0.07))
    K <- Kcom(cells, model=fit, restrict=TRUE)
    
    ##  Kscaled
    A <- Lscaled(japanesepines, renormalise=TRUE, correction="all")
  }
})
  
local({
  if(ALWAYS) {
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
    ## test validity without debugger,in case of quirks of compiler optimisation
    eZ <- edge.Ripley(Z, c(1,1))
    compere(eZ,      c(2,4/3), "at interior point of polygon (debug off)")
  }
})
#
# tests/kppm.R
#
# $Revision: 1.33 $ $Date: 2020/04/28 12:58:26 $
#
# Test functionality of kppm that depends on RandomFields
# Test update.kppm for old style kppm objects

local({

 fit <- kppm(redwood ~1, "Thomas") # sic
 fitx <- kppm(redwood ~x, "Thomas", verbose=TRUE)
 if(FULLTEST) {
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
   uu <- unitname(fitx)
   unitname(fitCx) <- "furlong"
   mo <- model.images(fitCx)
   p <- psib(fit)
   px <- psib(fitx)
 }
 if(ALWAYS) {
   Y <- simulate(fitx, seed=42, saveLambda=TRUE)[[1]]
 }

 if(FULLTEST) {
   #' vcov.kppm different algorithms
   vc  <- vcov(fitx)
   vc2 <- vcov(fitx, fast=TRUE)
   vc3 <- vcov(fitx, fast=TRUE, splitup=TRUE)
   vc4 <- vcov(fitx,            splitup=TRUE)

   ## other code blocks
   a <- varcount(fitx, function(x,y){x+1}) # always positive
   a <- varcount(fitx, function(x,y){y-1}) # always negative
   a <- varcount(fitx, function(x,y){x+y}) # positive or negative

   #' improve.kppm
   fitI <- update(fit, improve.type="quasi")
   fitxI <- update(fitx, improve.type="quasi")
   fitxIs <- update(fitx, improve.type="quasi", fast=FALSE) 
   #' vcov.kppm
   vcI <- vcov(fitxI)
 }

  ## plot.kppm including predict.kppm
 if(ALWAYS) {
   fitMC <- kppm(redwood ~ x, "Thomas")
   plot(fitMC)
 }
 if(FULLTEST) {
   fitCL <- kppm(redwood ~ x, "Thomas", method="c")
   fitPA <- kppm(redwood ~ x, "Thomas", method="p")
   plot(fitCL)
   plot(fitPA)

   ## fit with composite likelihood method [thanks to Abdollah Jalilian]
   fut <- kppm(redwood ~ x, "VarGamma", method="clik2", nu.ker=-3/8)
   kfut <- as.fv(fut)
 }
 
 if(require(RandomFields)) {
   fit0 <- kppm(redwood ~1, "LGCP")
   is.poisson(fit0)
   Y0 <- simulate(fit0, saveLambda=TRUE)[[1]]
   stopifnot(is.ppp(Y0))
   p0 <- psib(fit0) # issues a warning

   if(FULLTEST) {
     ## fit LGCP using K function: slow
     fit1 <- kppm(redwood ~x, "LGCP",
                  covmodel=list(model="matern", nu=0.3),
                  control=list(maxit=3))
     Y1 <- simulate(fit1, saveLambda=TRUE)[[1]]
     stopifnot(is.ppp(Y1))
   }

   ## fit LGCP using pcf
   fit1p <- kppm(redwood ~x, "LGCP",
                 covmodel=list(model="matern", nu=0.3),
                 statistic="pcf")
   Y1p <- simulate(fit1p, saveLambda=TRUE)[[1]]
   stopifnot(is.ppp(Y1p))

   ## .. and using different fitting methods
   if(FULLTEST) {
     fit1pClik <- update(fit1p, method="clik")
     fit1pPalm <- update(fit1p, method="palm")
   }
   ## image covariate (a different code block) 
   xx <- as.im(function(x,y) x, Window(redwood))
   fit1xx <- update(fit1p, . ~ xx, data=solist(xx=xx))
   Y1xx <- simulate(fit1xx, saveLambda=TRUE)[[1]]
   stopifnot(is.ppp(Y1xx))
   if(FULLTEST) {
     fit1xxVG <- update(fit1xx, clusters="VarGamma", nu=-1/4)
     Y1xxVG <- simulate(fit1xxVG, saveLambda=TRUE)[[1]]
     stopifnot(is.ppp(Y1xxVG))
   }
   fit1xxLG <- update(fit1xx, clusters="LGCP",
                      covmodel=list(model="matern", nu=0.3),
                      statistic="pcf")
   Y1xxLG <- simulate(fit1xxLG, saveLambda=TRUE, drop=TRUE)
   stopifnot(is.ppp(Y1xxLG))
   
   # ... and Abdollah's code
   if(FULLTEST) {
     fit2 <- kppm(redwood ~x, cluster="Cauchy", statistic="K")
     Y2 <- simulate(fit2, saveLambda=TRUE)[[1]]
     stopifnot(is.ppp(Y2))
   }
   # check package mechanism
   kraever("RandomFields")
   
 }

})

if(FULLTEST) {
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
}

if(FULLTEST) {
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
}

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

if(FULLTEST) {
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
}

reset.spatstat.options()
  
