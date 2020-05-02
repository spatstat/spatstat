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
#'   tests/tessera.R
#'   Tessellation code, not elsewhere tested
#'   $Revision: 1.8 $ $Date: 2020/05/02 01:32:58 $
#'
if(FULLTEST) {
local({
  W <- owin()
  Wsub <- square(0.5)
  X <- runifpoint(7, W)
  A <- dirichlet(X)
  marks(A) <- 1:nobjects(A)
  Z <- distmap(letterR, invert=TRUE)[letterR, drop=FALSE]
  H <- tess(xgrid=0:2, ygrid=0:3)
  #' discretisation of tiles
  V <- as.im(A)
  B <- tess(window=as.mask(W), tiles=tiles(A))
  #' logical images
  D <- tess(image=(Z > 0.2))
  U <- (Z > -0.2) # TRUE or NA
  E <- tess(image=U, keepempty=TRUE)
  G <- tess(image=U, keepempty=FALSE)
  #' methods
  flay <- function(op, ..., Rect=H, Poly=A, Img=E) {
    a <- do.call(op, list(Rect, ...))
    b <- do.call(op, list(Poly, ...))
    e <- do.call(op, list(Img, ...))
  }
  flay(reflect)
  flay(flipxy)
  flay(shift, vec=c(1,2))
  flay(scalardilate, f=2) 
  flay(rotate, angle=pi/3, centre=c(0, 0))
  flay(rotate, angle=pi/2)
  flay(affine, mat=matrix(c(1,2,0,1), 2, 2), vec=c(1,2))
  flay(affine, mat=diag(c(1,2)))
  flay(as.data.frame)
  ##
  unitname(A) <- "km"
  unitname(B) <- c("metre", "metres")
  unitname(B)
  print(B)
  Bsub <- B[c(3,5,7)]
  print(Bsub)
  tilenames(H) <- letters[seq_along(tilenames(H))]
  G <- tess(xgrid=(0:3)/3, ygrid=(0:3)/3)
  tilenames(G) <- letters[1:9]
  h <- tilenames(G)
  GG <- as.tess(tiles(G))
  #'
  Pe <- intersect.tess(A, Wsub, keepmarks=TRUE)
  Pm <- intersect.tess(A, as.mask(Wsub), keepmarks=TRUE)
  H <- dirichlet(runifpoint(4, W))
  AxH <- intersect.tess(A, H, keepmarks=TRUE) # A is marked, H is not
  HxA <- intersect.tess(H, A, keepmarks=TRUE) # A is marked, H is not
  
  b <- bdist.tiles(D)
  b <- bdist.tiles(A[c(3,5,7)])
  #'
  Eim <- as.im(E, W=letterR)
  #'
  #' chop.tess
  #'    horiz/vert lines
  W <- square(1)
  H <- infline(h=(2:4)/5)
  V <- infline(v=(3:4)/5)
  WH <- chop.tess(W, H)
  WV <- chop.tess(W, V)
  #'     polygonal tessellation
  D <- dirichlet(runifpoint(4))
  DH <- chop.tess(D, H)
  DV <- chop.tess(D, V)
  #'     image-based tessellation
  f <- function(x,y){factor(round(4* (x^2 + y^2)))}
  A <- tess(image=as.im(f, W=W))
  L <- infline(p=(1:3)/3, theta=pi/4)
  AL <- chop.tess(A, L)
  AH <- chop.tess(A, H)
  AV <- chop.tess(A, V)
  #'
  #' quantess
  #' quantess.owin
  a <- quantess(square(1), "x", 3)
  a <- quantess(square(1), "y", 3)
  a <- quantess(square(1), "rad", 5, origin=c(1/2, 1/3))
  a <- quantess(square(1), "ang", 7, origin=c(1/2, 1/3))
  ZFUN <- function(x,y){y-x}
  a <- quantess(square(1), ZFUN, 3)
  b <- quantess(letterR, "y", 3)
  #' quantess.ppp
  d <- quantess(cells, "y", 4)
  g <- quantess(demopat, "x", 5)
  g <- quantess(demopat, "y", 5)
  g <- quantess(demopat, "rad", 5, origin=c(4442, 4214))
  g <- quantess(demopat, "ang", 5, origin=c(4442, 4214))
  g <- quantess(demopat, ZFUN, 7)
  #' quantess.im
  D <- distmap(demopat)
  h <- quantess(D, "y", 4)
  h <- quantess(D, ZFUN, 5)
  g <- quantess(D, "rad", 5, origin=c(4442, 4214))
  g <- quantess(D, "ang", 5, origin=c(4442, 4214))
  #'
  X <- shift(chorley, vec = c(1e6, 0))
  tes <- quantess(X, "x", 4)
  if(anyDuplicated(tilenames(tes)))
    stop("quantess produced non-unique tilenames")
  ## 
  ## 
  da <- dirichletAreas(discretise(runifpoint(15, letterR)))
})
}
#
#   tests/testaddvar.R
#
# test addvar options
#
#   $Revision: 1.3 $  $Date: 2020/05/02 01:32:58 $

if(FULLTEST) {
local({
  X <-  rpoispp(function(x,y){exp(3+3*x)})
  model <- ppm(X ~y)
  addvar(model, "x", crosscheck=TRUE)
  addvar(model, "x", bw.input="quad")
  w <- square(0.5)
  addvar(model, "x", subregion=w)
  addvar(model, "x", subregion=w, bw.input="points")
  Z <- as.im(function(x,y) { x }, Window(X))
  addvar(model, Z)
})
}
#
#   tests/testparres.R
#
# additional test of parres
#
#  $Revision: 1.7 $  $Date: 2020/05/02 01:32:58 $
#

if(FULLTEST) {
local({
X <-  rpoispp(function(x,y){exp(3+x+2*x^2)})
model <- ppm(X ~x+y)

# options in parres (and code blocks in print.parres)
parres(model, "x")
parres(model, "x", smooth.effect=TRUE)
parres(model, "x", bw.input="quad")
w <- square(0.5)
parres(model, "x", subregion=w)
parres(model, "x", subregion=w, bw.input="quad")
f <- function(x,y) { x + y }
parres(model, f)

# check whether 'update.ppm' has messed up internals
mod2 <- update(model, ~x)
parres(mod2, "x")

#' other kinds of covariates
mod3 <- ppm(X ~ x + offset(y))
parres(mod3, "offset(y)")
Z <- distmap(runifpoint(3))
parres(mod3, Z)
mod4 <- ppm(X ~ sin(x), data=solist(B=Z))
parres(mod4, "sin(x)")
parres(mod4, "B")

#' models with interaction
mod5 <- ppm(cells ~ x, AreaInter(0.06))
parres(mod5, "x")
dlin <- distfun(copper$SouthLines)
copfit <- ppm(copper$SouthPoints ~ dlin, Geyer(1,1))
parres(copfit, "dlin")

#' covariate need not be specified if there is only one.
parres(mod5)
parres(copfit)

#' infrastructure
ltuae <- evalCovariate(42, cells)
LTUAE <- evalCovariate(ltuae, cells)

fit <- ppm(amacrine ~ x * marks, nd=16)
dmat <- model.depends(fit)
check.separable(dmat, "x", c(x=FALSE, marks=FALSE), FALSE)
check.separable(dmat, "x", c(FALSE, FALSE), FALSE)
check.separable(dmat, "x", c(x=FALSE, marks=TRUE), FALSE)
})
}
#'
#'     tests/threedee.R
#'
#'     Tests of 3D code 
#'
#'      $Revision: 1.8 $ $Date: 2020/05/02 01:32:58 $
#'

local({
  X <- runifpoint3(30)
  Y <- runifpoint3(20)
  if(FULLTEST) {
    A <- runifpoint3(10, nsim=2)
    Z <- ppsubset(X, 2:4)
  }
  ##
  if(ALWAYS) { # includes C code
    d <- pairdist(X, periodic=TRUE, squared=TRUE)
    d <- crossdist(X, Y, squared=TRUE)
    d <- crossdist(X, Y, squared=TRUE, periodic=TRUE)
    #' 
    h <- has.close(X, 0.2)
    h <- has.close(X, 0.2, periodic=TRUE)
    h <- has.close(X, 0.2, Y=Y)
    h <- has.close(X, 0.2, Y=Y, periodic=TRUE)
    #' code blocks not otherwise reached
    rmax <- 0.6 * max(nndist(X))
    g <- G3est(X, rmax=rmax, correction="rs")
    g <- G3est(X, rmax=rmax, correction="km")
    g <- G3est(X, rmax=rmax, correction="Hanisch")
    g <- G3est(X, rmax=rmax, sphere="ideal")
    g <- G3est(X, rmax=rmax, sphere="digital")
    v <- sphere.volume()
    v <- digital.volume()
    #' older code
    co <- coords(X)
    xx <- co$x
    yy <- co$y
    zz <- co$z
    gg1 <- g3engine(xx, yy, zz, correction="Hanisch G3")
    gg2 <- g3engine(xx, yy, zz, correction="minus sampling")
    ff1 <- f3engine(xx, yy, zz, correction="no")
    ff2 <- f3engine(xx, yy, zz, correction="minus sampling")
  }
  ##
  if(ALWAYS) {
    #'class support
    X <- runifpoint3(10)
    print(X)
    print(X %mark% runif(10))
    print(X %mark% factor(letters[c(1:5,5:1)]))
    print(X %mark% data.frame(a=1:10, b=runif(10)))
    da <- as.Date(paste0("2020-01-0", c(1:5,5:1)))
    print(X %mark% da)
    print(X %mark% data.frame(a=1:10, b=da))
  }
})
#'    tests/trigraph.R
#'
#'   Tests for C code in trigraf.c
#'   
#'  $Revision: 1.4 $  $Date: 2020/05/02 01:32:58 $
#'
if(ALWAYS) { # depends on C code 
local({
  #' called from deldir.R
  spatstat.deldir.setopt(FALSE, TRUE)
  A <- delaunay(redwood)
  spatstat.deldir.setopt(FALSE, FALSE)
  B <- delaunay(redwood)
  spatstat.deldir.setopt(TRUE, TRUE)
  #' called from edges2triangles.R
  tryangles <- function(iedge, jedge, nt=0) {
    spatstat.options(fast.trigraph=FALSE)
    A <- edges2triangles(iedge, jedge)
    spatstat.options(fast.trigraph=TRUE)
    B <- edges2triangles(iedge, jedge)
    if(!all(dim(A) == dim(B)) || !all(A == B))
      stop(paste("Discrepancy in edges2triangles (with", nt, "triangles)"))
  }
  ii <- simplenet$from
  jj <- simplenet$to
  tryangles(ii,          jj,          0)
  tryangles(c(ii, 1),    c(jj, 5),    1)
  tryangles(c(ii, 1, 8), c(jj, 5, 9), 2)
})
}
reset.spatstat.options()


#
# tests/triplets.R
#
# test code for triplet interaction and associated summary function Tstat 
#
# $Revision: 1.8 $ $Date: 2020/05/02 01:32:58 $
#
if(ALWAYS) { # C code, platform dependence
local({
  #' valid model
  fit <- ppm(cells ~1, Triplets(0.1))
  fit
  suffstat(fit)
  #' invalid model 
  fitR <- ppm(redwood ~1, Triplets(0.1))
  fitR
  suffstat(fitR)
  #' hard core (zero triangles, coefficient is NA)
  fit0 <- ppm(cells ~1, Triplets(0.05))
  fit0
  suffstat(fit0)
  #' bug case (1 triangle in data)
  fit1 <- ppm(cells ~1, Triplets(0.15))
  fit1
  suffstat(fit1)
  #' Tstat function, all code blocks
  a <- Tstat(redwood, ratio=TRUE,
             correction=c("none", "border", "bord.modif", "translate"))
  #' simulation
  X <- simulate(fit)
  mod <- list(cif="triplets",par=list(beta=50,gamma=0.2,r=0.07), w=square(1))
  Xm <- rmh(model=mod,start=list(n.start=5), control=list(nrep=1e5))
  #' hard core
  mod$par$gamma <- 0
  XmHard <- rmh(model=mod,start=list(n.start=5), control=list(nrep=1e5))
})
}
