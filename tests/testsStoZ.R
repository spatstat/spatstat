#
#  tests/segments.R
#
#  $Revision: 1.11 $  $Date: 2017/02/20 10:15:30 $

require(spatstat)

local({
# pointed out by Jeff Laake
W <- owin()
X <- psp(x0=.25,x1=.25,y0=0,y1=1,window=W)
X[W]

# migrated from 'lpp'

X <- psp(runif(10),runif(10),runif(10),runif(10), window=owin())
Z <- as.mask.psp(X)
Z <- pixellate(X)

# more tests of lppm code

fit <- lppm(unmark(chicago) ~ polynom(x,y,2))
Z <- predict(fit)

# tests of pixellate.psp -> seg2pixL

ns <- 50
out <- numeric(ns)
for(i in 1:ns) {
  X <- psp(runif(1), runif(1), runif(1), runif(1), window=owin())
  len <- lengths.psp(X)
  dlen <- sum(pixellate(X)$v)
  out[i] <- if(len > 1e-7) dlen/len else 1
}
if(diff(range(out)) > 0.01) stop(paste(
       "pixellate.psp test 1: relative error [",
       paste(diff(range(out)), collapse=", "),
       "]"))

# Michael Sumner's test examples

set.seed(33)
n <- 2001
co <- cbind(runif(n), runif(n))
ow <- owin()
X <- psp(co[-n,1], co[-n,2], co[-1,1], co[-1,2], window=ow)
s1 <- sum(pixellate(X))
s2 <- sum(lengths.psp(X))
if(abs(s1 - s2)/s2 > 0.01) {
  stop(paste("pixellate.psp test 2:",
             "sum(pixellate(X)) = ", s1,
             "!=", s2, "= sum(lengths.psp(X))"))
}

wts <- 1/(lengths.psp(X) * X$n)
s1 <- sum(pixellate(X, weights=wts))
if(abs(s1-1) > 0.01) {
  stop(paste("pixellate.psp test 3:",
             "sum(pixellate(X, weights))=", s1,
             " (should be 1)"))
}

X <- psp(0, 0, 0.01, 0.001, window=owin())
s1 <- sum(pixellate(X))
s2 <- sum(lengths.psp(X))
if(abs(s1 - s2)/s2 > 0.01) {
  stop(paste("pixellate.psp test 4:",
             "sum(pixellate(X)) = ", s1,
             "!=", s2, "= sum(lengths.psp(X))"))
}

X <- psp(0, 0, 0.001, 0.001, window=owin())
s1 <- sum(pixellate(X))
s2 <- sum(lengths.psp(X))
if(abs(s1 - s2)/s2 > 0.01) {
  stop(paste("pixellate.psp test 5:",
             "sum(pixellate(X)) = ", s1,
             "!=", s2, "= sum(lengths.psp(X))"))
}

#' tests of density.psp
Y <- as.psp(simplenet)
YC <- density(Y, 0.2, method="C", edge=FALSE, dimyx=64)
YI <- density(Y, 0.2, method="interpreted", edge=FALSE, dimyx=64)
YF <- density(Y, 0.2, method="FFT", edge=FALSE, dimyx=64)
xCI <- max(abs(YC/YI - 1))
xFI <- max(abs(YF/YI - 1))
if(xCI > 0.01) stop(paste("density.psp C algorithm relative error =", xCI))
if(xFI > 0.01) stop(paste("density.psp FFT algorithm relative error =", xFI))

})
#
## tests/sigtraceprogress.R
#
## Tests of *.sigtrace and *.progress
#
## $Revision: 1.3 $ $Date: 2015/12/29 08:54:49 $

require(spatstat)
local({
  plot(dclf.sigtrace(redwood, nsim=19, alternative="greater", rmin=0.02,
                     verbose=FALSE))
  plot(dclf.progress(redwood, nsim=19, alternative="greater", rmin=0.02,
                     verbose=FALSE))
  plot(dg.sigtrace(redwood, nsim=5, alternative="greater", rmin=0.02,
                     verbose=FALSE))
  plot(dg.progress(redwood, nsim=5, alternative="greater", rmin=0.02,
                   verbose=FALSE))
  ## test 'leave-two-out' algorithm
  a <- dclf.sigtrace(redwood, Lest, nsim=9, use.theory=FALSE, leaveout=2,
                     verbose=FALSE)
  aa <- dclf.progress(redwood, Lest, nsim=9, use.theory=FALSE, leaveout=2,
                      verbose=FALSE)
  b <- dg.sigtrace(redwood, Lest, nsim=5, use.theory=FALSE, leaveout=2)
  bb <- dg.progress(redwood, Lest, nsim=5, use.theory=FALSE, leaveout=2,
                    verbose=FALSE)
})
#
# tests/slrm.R
#
# $Revision: 1.1 $ $Date: 2013/04/19 10:14:52 $
#
# Test slrm fitting and prediction when there are NA's
#

require(spatstat)
local({
  X <- copper$SouthPoints
  W <- owin(poly=list(x=c(0,35,35,1),y=c(1,1,150,150)))
  Y <- X[W]
  fit <- slrm(Y ~ x+y)
  pred <- predict(fit)
})


#'    tests/sparse3Darrays.R
#'  Basic tests of sparse3array.R code
#'  $Revision: 1.8 $ $Date: 2017/02/22 09:00:27 $

require(spatstat)
local({

  if(require(Matrix)) {
    M <- sparse3Darray(i=1:4, j=sample(1:4, replace=TRUE),
                       k=c(1,2,1,2), x=1:4, dims=c(5,5,2))

    M

    dimnames(M) <- list(letters[1:5], LETTERS[1:5], c("yes", "no"))
    M
    
    U <- aperm(M, c(1,3,2))
    U
    
    M[ 3:4, , ]
    
    M[ 3:4, 2:4, ]
    
    M[, 3, ]

    M[, 3, , drop=FALSE]
    
    MA <- as.array(M)
    UA <- as.array(U)

    ## tests of "[<-.sparse3Darray"
    Mflip <- Mzero <- MandM <- M
    Mflip[ , , 2:1] <- M
    stopifnot(Mflip[3,1,1] == M[3,1,2])
    Mzero[1:3,1:3,] <- 0
    stopifnot(all(Mzero[1,1,] == 0))
    M2a <- M[,,2,drop=FALSE]
    M2d <- M[,,2,drop=TRUE]
    MandM[,,1] <- M2a
    MandM[,,1] <- M2d

    # matrix index
    M[cbind(3:5, 2, 2)]
    
    ## tests of arithmetic (Math, Ops, Summary)
    negM <- -M
    oneM <- 1 * M
    twoM <- M + M
    range(M)

    cosM <- cos(M)  # non-sparse
    sinM <- sin(M)  # sparse
    
    stopifnot(all((M+M) == 2*M))     # non-sparse
    stopifnot(!any((M+M) != 2*M))    # sparse

    ztimesM <- (1:5) * M  # sparse
    zplusM <- (1:5) + M  # non-sparse
    
    ## tensor operator

    tenseur(c(1,-1), M, 1, 3)
    tenseur(M, M, 1:2, 1:2)
    tenseur(M, M, 1:2, 2:1)
    V <- sparseVector(i=c(1,3,6),x=1:3, length=7)
    tenseur(V,V)
    tenseur(V,V,1,1)

    ## test of anyNA method
    anyNA(M)
    
    ## a possible application in spatstat
    cl10 <- as.data.frame(closepairs(cells, 0.1))
    cl12 <- as.data.frame(closepairs(cells, 0.12))
    cl10$k <- 1
    cl12$k <- 2
    cl <- rbind(cl10, cl12)
    n <- npoints(cells)
    Z <- with(cl,
              sparse3Darray(i=i, j=j, k=k, x=1, dims=c(n,n,2)))
    dimnames(Z) <- list(NULL, NULL, c("r=0.1", "r=0.12"))

    Z <- aperm(Z, c(3,1,2))
    stopifnot(all(sumsymouterSparse(Z) == sumsymouter(as.array(Z))))

    # no entries indexed
    Z[integer(0), integer(0), integer(0)] <- 42
    Z[matrix(, 0, 3)] <- 42
  }
})


#
#  tests/splitpea.R
#
#  Check behaviour of split.ppp etc
#
#  Thanks to Marcelino de la Cruz
#
#  $Revision: 1.11 $  $Date: 2016/03/05 01:33:47 $
#

require(spatstat)

local({
W <- square(8)
X <- ppp(c(2.98, 4.58, 7.27, 1.61, 7.19),
         c(7.56, 5.29, 5.03, 0.49, 1.65),
         window=W)
Z <- quadrats(W, 4, 4)
Yall <- split(X, Z, drop=FALSE)
Ydrop <- split(X, Z, drop=TRUE)

P <- Yall[[1]]
if(!all(inside.owin(P$x, P$y, P$window)))
  stop("Black hole detected when drop=FALSE")
P <- Ydrop[[1]]
if(!all(inside.owin(P$x, P$y, P$window)))
  stop("Black hole detected when drop=TRUE")

Ydrop[[1]] <- P[1]
split(X, Z, drop=TRUE) <- Ydrop

# test NA handling
Zbad <- quadrats(square(4), 2, 2)
Ybdrop <- split(X, Zbad, drop=TRUE)
Yball  <- split(X, Zbad, drop=FALSE)

# From Marcelino
set.seed(1)
W<- square(10) # the big window
puntos<- rpoispp(0.5, win=W)
data(letterR)
r00 <- letterR
r05 <- shift(letterR,c(0,5))
r50 <- shift(letterR,c(5,0))
r55 <- shift(letterR,c(5,5))
tessr4 <- tess(tiles=list(r00, r05,r50,r55))
puntosr4 <- split(puntos, tessr4, drop=TRUE)
split(puntos, tessr4, drop=TRUE) <- puntosr4

## More headaches with mark format
A <- runifpoint(10)
B <- runifpoint(10)
AB <- split(superimpose(A=A, B=B))

#' check that split<- respects ordering where possible
X <- amacrine
Y <- split(X)
split(X) <- Y
stopifnot(identical(X, amacrine))

#' split.ppx
df <- data.frame(x=runif(4),y=runif(4),t=runif(4),
                 age=rep(c("old", "new"), 2),
                 mineral=factor(rep(c("Au","Cu"), each=2),
                                levels=c("Au", "Cu", "Pb")),
                 size=runif(4))
X <- ppx(data=df, coord.type=c("s","s","t","m", "m","m"))
Y <- split(X, "age")
Y <- split(X, "mineral", drop=TRUE)

})
#
#   tests/step.R
#
#   $Revision: 1.4 $  $Date: 2015/12/29 08:54:49 $
#
# test for step() operation
#
require(spatstat)
local({
  Z <- as.im(function(x,y){ x^3 - y^2 }, nztrees$window)
  fitP <- ppm(nztrees ~x+y+Z, covariates=list(Z=Z))
  step(fitP)
  fitS <- update(fitP, Strauss(7))
  step(fitS)
  fitM <- ppm(amacrine ~ marks*(x+y),
              MultiStrauss(types=levels(marks(amacrine)), radii=matrix(0.04, 2, 2)))
  step(fitM)
})


##
## tests/symbolmaps.R
##
##   Quirks associated with symbolmaps, etc.
##
## $Revision: 1.3 $ $Date: 2015/12/29 08:54:49 $

local({
  require(spatstat)
  set.seed(100)
  
  ## spacing too large for tiles - upsets various pieces of code
  V <- as.im(dirichlet(runifpoint(8)))
  textureplot(V, spacing=2)

  g1 <- symbolmap(range=c(0,100), size=function(x) x/50)
  invoke.symbolmap(g1, 50, x=numeric(0), y=numeric(0), add=TRUE)

})
#
#   tests/testaddvar.R
#
# test addvar options
#
#   $Revision: 1.2 $  $Date: 2015/12/29 08:54:49 $

X <-  rpoispp(function(x,y){exp(3+3*x)})
model <- ppm(X ~y)
addvar(model, "x", crosscheck=TRUE)
addvar(model, "x", bw.input="quad")
w <- square(0.5)
addvar(model, "x", subregion=w)
addvar(model, "x", subregion=w, bw.input="points")
#
#   tests/testparres.R
#
# additional test of parres
#
#  $Revision: 1.2 $  $Date: 2015/12/29 08:54:49 $
#
require(spatstat)
local({
X <-  rpoispp(function(x,y){exp(3+x+2*x^2)})
model <- ppm(X ~x+y)

# options in parres
parres(model, "x")
parres(model, "x", bw.input="quad")
w <- square(0.5)
parres(model, "x", subregion=w)
parres(model, "x", subregion=w, bw.input="quad")

# check whether 'update.ppm' has messed up internals
mod2 <- update(model, ~x)
parres(mod2, "x")
})
#
# tests/triplets.R
#
# test code for triplet interaction
#
# $Revision: 1.5 $ $Date: 2015/12/29 08:54:49 $
#
require(spatstat)
local({
  fit <- ppm(redwood ~1, Triplets(0.1))
  fit
  suffstat(fit)
  # hard core (zero triangles, coefficient is NA)
  fit0 <- ppm(cells ~1, Triplets(0.05))
  fit0
  suffstat(fit0)
  # bug case (1 triangle in data)
  fit1 <- ppm(cells ~1, Triplets(0.15))
  fit1
  suffstat(fit1)
})
#
#  tests/undoc.R
#
#   $Revision: 1.3 $   $Date: 2017/02/20 10:51:56 $
#
#  Test undocumented hacks, etc

require(spatstat)
local({
  # pixellate.ppp accepts a data frame of weights
  pixellate(cells, weights=data.frame(a=1:42, b=42:1))
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
#  $Revision: 1.6 $  $Date: 2015/12/29 08:54:49 $
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
  modHyb <- ppm(japanesepines ~ 1, Hybrid(Strauss(0.05), Strauss(0.1)))
})
#
# tests/windows.R
#
# Tests of owin geometry code
#
#  $Revision: 1.3 $  $Date: 2015/12/29 08:54:49 $

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
})

##
## tests/xysegment.R
##
##    Test weird problems and boundary cases for line segment code
##
##    $Version$ $Date: 2016/02/12 08:18:08 $ 
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
})
