#
#  tests/undoc.R
#
#   $Revision: 1.11 $   $Date: 2019/11/29 03:41:54 $
#
#  Test undocumented hacks, experimental code, etc

require(spatstat)
local({
  ## cases of 'pickoption'
  aliases <- c(Lenin="Ulyanov", Stalin="Djugashvili", Trotsky="Bronstein")
  surname <- "Trot"
  pickoption("leader", surname,  aliases)
  pickoption("leader", surname,  aliases, exact=TRUE, die=FALSE)
  ## pixellate.ppp accepts a data frame of weights
  pixellate(cells, weights=data.frame(a=1:42, b=42:1))
  ## test parts of 'rmhsnoop' that don't require interaction with user
  rmhSnoopEnv(cells, Window(cells), 0.1)
  ## linim helper functions
  df <- pointsAlongNetwork(simplenet, 0.2)
  ## Berman-Turner frame
  A <- bt.frame(quadscheme(cells), ~x, Strauss(0.07), rbord=0.07)
  print(A)
  ## digestCovariates
  D <- distfun(cells)
  Z <- distmap(cells)
  U <- dirichlet(cells)
  stopifnot(is.scov(D))
  stopifnot(is.scov(Z))
  stopifnot(is.scov(U))
  stopifnot(is.scov("x"))
  dg <- digestCovariates(D=D,Z=Z,U=U,"x",list(A="x", B=D))
  ##
  a <- getfields(dg, c("A", "D", "niets"), fatal=FALSE)
  ## util.R
  gg <- pointgrid(owin(), 7)
  checkbigmatrix(1000000L, 1000000L, FALSE, TRUE)
  spatstatDiagnostic("whatever")
  M <- list(list(a=2, b=FALSE),
            list(a=2, b=TRUE))
  stopifnot(!allElementsIdentical(M))
  stopifnot(allElementsIdentical(M, "a"))
  ##
  A <- Strauss(0.1)
  A <- reincarnate.interact(A)
  ##
  ## special lists
  B <- solist(a=cells, b=redwood, c=japanesepines)
  BB <- as.ppplist(B)
  BL <- as.layered(B)
  DB <- as.imlist(lapply(B, density))
  is.solist(B)
  is.ppplist(B)
  is.imlist(DB)
  ## case of density.ppplist 
  DEB <- density(BB, se=TRUE)
  
  ## fft
  z <- matrix(1:16, 4, 4)
  a <- fft2D(z, west=FALSE)
  if(fftwAvailable())
    b <- fft2D(z, west=TRUE)

  ## experimental interactions
  pot <- function(d, par) { d <= 0.1 }
  A <- Saturated(pot)
  print(A)
  A <- update(A, name="something")
  ppm(amacrine ~ x, A, rbord=0.1)

  ## infrastructure of iplot, istat
  fakepanel <- list(x=redwood,
                    xname="redwood",
                    envel="none",
                    stat="data",
                    sigma=0.08,
                    pcfbw=0.01,
                    simx=rpoispp(ex=redwood, nsim=39))
  a <- do.istat(fakepanel)
  for(envel in c("none", "pointwise", "simultaneous")) {
    fakepanel$envel <- envel
    for(stat in c("density",
                  "Kest", "Lest", "pcf",
                  "Kinhom", "Linhom", "Fest", "Gest", "Jest")) {
      fakepanel$stat <- stat
      a <- do.istat(fakepanel)
    }
  }

  #' version-checking
  now <- Sys.Date()
  versioncurrency.spatstat(now + 80, FALSE)
  versioncurrency.spatstat(now + 140, FALSE)
  versioncurrency.spatstat(now + 400, FALSE)
  versioncurrency.spatstat(now + 1000)
})

##
##  tests/updateppm.R
##
##  Check validity of update.ppm
##
##  $Revision: 1.4 $ $Date: 2016/03/08 06:30:46 $

local({
    require(spatstat)
    h <- function(m1, m2) {
        mc <- deparse(sys.call())
        cat(paste(mc, "\t... "))
        m1name <- deparse(substitute(m1))
        m2name <- deparse(substitute(m2))
        if(!identical(names(coef(m1)), names(coef(m2))))
            stop(paste("Differing results for", m1name, "and", m2name,
                       "in updateppm.R"),
                 call.=FALSE)
        cat("OK\n")
    }
    X <- redwood[c(TRUE,FALSE)]
    Y <- redwood[c(FALSE,TRUE)]
    fit0f <- ppm(X ~ 1, nd=8)
    fit0p <- ppm(X, ~1, nd=8)
    fitxf <- ppm(X ~ x, nd=8)
    fitxp <- ppm(X, ~x, nd=8)

    cat("Basic consistency ...\n")
    h(fit0f, fit0p)
    h(fitxf, fitxp)

    cat("\nTest correct handling of model formulas ...\n")
    h(update(fitxf, Y), fitxf)
    h(update(fitxf, Q=Y), fitxf)
    h(update(fitxf, Y~x), fitxf)
    h(update(fitxf, Q=Y~x), fitxf)
    h(update(fitxf, ~x), fitxf)

    h(update(fitxf, Y~1), fit0f)
    h(update(fitxf, ~1), fit0f)
    h(update(fit0f, Y~x), fitxf)
    h(update(fit0f, ~x), fitxf)

    h(update(fitxp, Y), fitxp)
    h(update(fitxp, Q=Y), fitxp)
    h(update(fitxp, Y~x), fitxp)
    h(update(fitxp, Q=Y~x), fitxp)
    h(update(fitxp, ~x), fitxp)

    h(update(fitxp, Y~1), fit0p)
    h(update(fitxp, ~1), fit0p)
    h(update(fit0p, Y~x), fitxp)
    h(update(fit0p, ~x), fitxp)

    cat("\nTest scope handling for left hand side ...\n")
    X <- Y
    h(update(fitxf), fitxf)

    cat("\nTest scope handling for right hand side ...\n")
    Z <- distmap(X)
    fitZf <- ppm(X ~ Z)
    fitZp <- ppm(X, ~ Z)
    h(update(fitxf, X ~ Z), fitZf)
    h(update(fitxp, X ~ Z), fitZp)
    h(update(fitxf, . ~ Z), fitZf)
    h(update(fitZf, . ~ x), fitxf)
    h(update(fitZf, . ~ . - Z), fit0f)
    h(update(fitxp, . ~ Z), fitZp)
    h(update(fitZp, . ~ . - Z), fit0p)
    h(update(fit0p, . ~ . + Z), fitZp)
    h(update(fitZf, . ~ . ), fitZf)
    h(update(fitZp, . ~ . ), fitZp)

    cat("\nTest use of internal data ...\n")
    h(update(fitZf, ~ x, use.internal=TRUE), fitxf)
    fitsin <- update(fitZf, X~sin(Z))
    h(update(fitZf, ~ sin(Z), use.internal=TRUE), fitsin)

    cat("\nTest step() ... ")
    fut <- ppm(X ~ Z + x + y, nd=8)
    fut0 <- step(fut, trace=0)
    cat("OK\n")
})

# test update.lppm

local({
  X <- runiflpp(20, simplenet)
  fit0 <- lppm(X ~ 1)
  fit1 <- update(fit0, ~ x)
  anova(fit0, fit1, test="LR")
  cat("update.lppm(fit, ~trend) is OK\n")
  fit2 <- update(fit0, . ~ x)
  anova(fit0, fit2, test="LR")
  cat("update.lppm(fit, . ~ trend) is OK\n")
})
#
#  tests/vcovppm.R
#
#  Check validity of vcov.ppm algorithms
#
#  Thanks to Ege Rubak
#
#  $Revision: 1.11 $  $Date: 2020/01/11 09:38:40 $
#

require(spatstat)

local({

  set.seed(42)
  X <- rStrauss(200, .5, .05)
  model <- ppm(X, inter = Strauss(.05))

  b  <- vcov(model, generic = TRUE, algorithm = "basic")
  v  <- vcov(model, generic = TRUE, algorithm = "vector")
  vc <- vcov(model, generic = TRUE, algorithm = "vectorclip")
  vn <- vcov(model, generic = FALSE)

  disagree <- function(x, y, tol=1e-7) { max(abs(x-y)) > tol }
  asymmetric <- function(x) { disagree(x, t(x)) }

  if(asymmetric(b))
    stop("Non-symmetric matrix produced by vcov.ppm 'basic' algorithm")
  if(asymmetric(v))
    stop("Non-symmetric matrix produced by vcov.ppm 'vector' algorithm")
  if(asymmetric(vc))
    stop("Non-symmetric matrix produced by vcov.ppm 'vectorclip' algorithm")
  if(asymmetric(vn))
    stop("Non-symmetric matrix produced by vcov.ppm Strauss algorithm")
    
  if(disagree(v, b))
    stop("Disagreement between vcov.ppm algorithms 'vector' and 'basic' ")
  if(disagree(v, vc))
    stop("Disagreement between vcov.ppm algorithms 'vector' and 'vectorclip' ")
  if(disagree(vn, vc))
    stop("Disagreement between vcov.ppm generic and Strauss algorithms")

  # Geyer code
  xx <- c(0.7375956, 0.6851697, 0.6399788, 0.6188382)
  yy <- c(0.5816040, 0.6456319, 0.5150633, 0.6191592)
  Y <- ppp(xx, yy, window=square(1))
  modelY <- ppm(Y ~1, Geyer(0.1, 1))

  b  <- vcov(modelY, generic = TRUE, algorithm = "basic")
  v  <- vcov(modelY, generic = TRUE, algorithm = "vector")
  vc <- vcov(modelY, generic = TRUE, algorithm = "vectorclip")

  if(asymmetric(b))
    stop("Non-symmetric matrix produced by vcov.ppm 'basic' algorithm for Geyer model")
  if(asymmetric(v))
    stop("Non-symmetric matrix produced by vcov.ppm 'vector' algorithm for Geyer model")
  if(asymmetric(vc))
    stop("Non-symmetric matrix produced by vcov.ppm 'vectorclip' algorithm for Geyer model")
  
  if(disagree(v, b))
    stop("Disagreement between vcov.ppm algorithms 'vector' and 'basic' for Geyer model")
  if(disagree(v, vc))
    stop("Disagreement between vcov.ppm algorithms 'vector' and 'vectorclip' for Geyer model")

  ## tests of 'deltasuffstat' code
  ##     Handling of offset terms
  modelH <- ppm(cells ~x, Hardcore(0.05))
  a <- vcov(modelH, generic=TRUE) ## may fall over
  b <- vcov(modelH, generic=FALSE)
  if(disagree(a, b))
    stop("Disagreement between vcov.ppm algorithms for Hardcore model")
  
  ##     Correctness of pairwise.family$delta2
  modelZ <- ppm(amacrine ~1, MultiStrauss(radii=matrix(0.1, 2, 2)))
  b <- vcov(modelZ, generic=FALSE)
  g <- vcov(modelZ, generic=TRUE)
  if(disagree(b, g))
    stop("Disagreement between vcov.ppm algorithms for MultiStrauss model")

  ## Test that 'deltasuffstat' works for Hybrids
  modelHyb <- ppm(japanesepines ~ 1, Hybrid(Strauss(0.05), Strauss(0.1)))
  vHyb <- vcov(modelHyb)

  ## Code blocks for other choices of 'what'
  model <- ppm(X ~ 1, Strauss(.05))
  cG <- vcov(model, what="corr")
  cP <- vcov(update(model, Poisson()), what="corr")
  ## outdated usage
  cX <- vcov(model, A1dummy=TRUE)

  ## Model with zero-length coefficient vector
  lam <- intensity(X)
  f <- function(x,y) { rep(lam, length(x)) }
  model0 <- ppm(X ~ offset(log(f)) - 1)
  dd <- vcov(model0)
  cc <- vcov(model0, what="corr")

  ## Model with NA coefficients
  fit <- ppm(X ~ log(f))
  vcov(fit)
  fitE <- emend(fit, trace=TRUE)
  
  ## Other weird stuff
  su <- suffloc(ppm(X ~ x))
})
#
# tests/windows.R
#
# Tests of owin geometry code
#
#  $Revision: 1.14 $  $Date: 2020/01/23 04:02:20 $

require(spatstat)
local({
  # Ege Rubak spotted this problem in 1.28-1
  A <- as.owin(ants)
  B <- dilation(A, 140)
  if(!is.subset.owin(A, B))
    stop("is.subset.owin fails in polygonal case")

  # thanks to Tom Rosenbaum
  A <- shift(square(3), origin="midpoint")
  B <- shift(square(1), origin="midpoint")
  AB <- setminus.owin(A, B)
  D <- shift(square(2), origin="midpoint")
  if(is.subset.owin(D,AB))
    stop("is.subset.owin fails for polygons with holes")

  ## thanks to Brian Ripley / SpatialVx
  M <- as.mask(letterR)
  stopifnot(area(bdry.mask(M)) > 0)
  stopifnot(area(convexhull(M)) > 0)
  R <- as.mask(square(1))
  stopifnot(area(bdry.mask(R)) > 0)
  stopifnot(area(convexhull(R)) > 0)

  ##
  RR <- convexify(as.mask(letterR))
  CC <- covering(letterR, 0.05, eps=0.1)
  
  #' as.owin.data.frame
  V <- as.mask(letterR, eps=0.2)
  Vdf <- as.data.frame(V)
  Vnew <- as.owin(Vdf)
  zz <- mask2df(V)

  RM <- owinpoly2mask(letterR, as.mask(Frame(letterR)), check=TRUE)

  #' as.owin
  U <- as.owin(quadscheme(cells))
  U2 <- as.owin(list(xmin=0, xmax=1, ymin=0, ymax=1))
  
  #' intersections involving masks
  B1 <- square(1)
  B2 <- as.mask(shift(B1, c(0.2, 0.3)))
  o12 <- overlap.owin(B1, B2)
  o21 <- overlap.owin(B2, B1)
  i12 <- intersect.owin(B1, B2, eps=0.01)
  i21 <- intersect.owin(B2, B1, eps=0.01)
  E2 <- emptywindow(square(2))
  e12 <- intersect.owin(B1, E2)
  e21 <- intersect.owin(E2, B1)
  
  #' geometry
  inradius(B1)
  inradius(B2)
  inradius(letterR)
  inpoint(B1)
  inpoint(B2)
  inpoint(letterR)
  is.convex(B1)
  is.convex(B2)
  is.convex(letterR)
  volume(letterR)
  perimeter(as.mask(letterR))
  boundingradius(cells)
  
  boundingbox(letterR)
  boundingbox(letterR, NULL)
  boundingbox(cells, ppm(cells ~ 1))
  boundingbox(solist(letterR))

  spatstat.options(Cbdrymask=FALSE)
  bb <- bdry.mask(letterR)
  spatstat.options(Cbdrymask=TRUE)
  
  X <- longleaf[square(50)]
  marks(X) <- marks(X)/8
  D <- discs(X, delta=5, separate=TRUE)

  AD <- dilated.areas(cells,
                      r=0.01 * matrix(1:10, 10,1),
                      constrained=FALSE, exact=FALSE)
  
  periodify(B1, 2)
  periodify(union.owin(B1, B2), 2)
  periodify(letterR, 2)

  #' Ancient bug in inside.owin
  W5 <- owin(poly=1e5*cbind(c(-1,1,1,-1),c(-1,-1,1,1)))
  W6 <- owin(poly=1e6*cbind(c(-1,1,1,-1),c(-1,-1,1,1)))
  i5 <- inside.owin(0,0,W5)
  i6 <- inside.owin(0,0,W6)
  if(!i5) stop("Wrong answer from inside.owin")
  if(i5 != i6) stop("Results from inside.owin are scale-dependent")
  
  #' miscellaneous utilities
  thrash <- function(f) {
    f(letterR)
    f(Frame(letterR))
    f(as.mask(letterR))
  }
  thrash(meanX.owin)
  thrash(meanY.owin)
  thrash(intX.owin)
  thrash(intY.owin)

  interpretAsOrigin("right", letterR)
  interpretAsOrigin("bottom", letterR)
  interpretAsOrigin("bottomright", letterR)
  interpretAsOrigin("topleft", letterR)
  interpretAsOrigin("topright", letterR)

  A <- break.holes(letterR)
  B <- break.holes(letterR, splitby="y")
  plot(letterR, col="blue", use.polypath=FALSE)
})


local({
  #' mask conversion
  M <- as.mask(letterR)
  D2 <- as.data.frame(M)              # two-column
  D3 <- as.data.frame(M, drop=FALSE)  # three-column
  M2 <- as.owin(D2)
  M3 <- as.owin(D3)
  W2 <- owin(mask=D2)
  W3 <- owin(mask=D3)
  #' void/empty cases
  nix <- nearest.raster.point(numeric(0), numeric(0), M)
  E <- emptywindow(Frame(letterR))
  print(E)
  #' cases of summary.owin
  print(summary(E)) # empty
  print(summary(Window(humberside))) # single polygon
  #' additional cases of owin()
  B <- owin(mask=M$m) # no pixel size or coordinate info
  xy <- as.data.frame(letterR)
  xxyy <- split(xy[,1:2], xy$id)
  spatstat.options(checkpolygons=TRUE)
  H <- owin(poly=xxyy, check=TRUE)
})

reset.spatstat.options()
##
## tests/xysegment.R
##                      [SEE ALSO tests/segments.R]
##
##    Test weird problems and boundary cases for line segment code
##
##    $Version$ $Date: 2018/05/13 04:22:28 $ 
##
require(spatstat)
local({
  # segment of length zero
  B <- psp(1/2, 1/2, 1/2, 1/2, window=square(1))
  BB <- angles.psp(B)
  A <- runifpoint(3)
  AB <- project2segment(A,B)

  # mark inheritance
  X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  marks(X) <- 1:10
  Y <- selfcut.psp(X)
  marks(X) <- data.frame(A=1:10, B=factor(letters[1:10]))
  Z <- selfcut.psp(X)

  #' psp class support
  S <- as.psp(simplenet)
  marks(S) <- sample(factor(c("A","B")), nobjects(S), replace=TRUE)
  intensity(S)
  intensity(S, weights=runif(nsegments(S)))
})
