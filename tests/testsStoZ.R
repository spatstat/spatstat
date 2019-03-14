#'   tests/sdr.R
#'
#'   $Revision: 1.1 $ $Date: 2018/05/13 03:14:49 $

require(spatstat)
local({
  AN <- sdr(bei, bei.extra, method="NNIR")
  AV <- sdr(bei, bei.extra, method="SAVE")
  AI <- sdr(bei, bei.extra, method="SIR")
  AT <- sdr(bei, bei.extra, method="TSE")
  subspaceDistance(AN$B, AV$B)
  dimhat(AN$M)
})
#
#  tests/segments.R
#
#  $Revision: 1.13 $  $Date: 2018/07/22 02:15:02 $

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

#' misc
PX <- periodify(X, 2)

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

#' as.psp.data.frame

  df <- as.data.frame(matrix(runif(40), ncol=4))
  A <- as.psp(df, window=square(1))
  colnames(df) <- c("x0","y0","x1","y1")
  df <- cbind(df, data.frame(marks=1:nrow(df)))
  B <- as.psp(df, window=square(1))
  colnames(df) <- c("xmid", "ymid", "length", "angle", "marks")
  E <- as.psp(df, window=square(c(-1,2)))

#' print and summary methods
  A
  B
  E
  summary(B)
  M <- B
  marks(M) <- data.frame(id=marks(B), len=lengths.psp(B))
  M
  summary(M)

#' undocumented  
  as.ppp(B)

  #' segment crossing code
  X <- psp(runif(30),runif(30),runif(30),runif(30), window=owin())
  spatstat.options(selfcrossing.psp.useCall=FALSE)
  selfcrossing.psp(X)
  spatstat.options(selfcrossing.psp.useCall=TRUE)
})

reset.spatstat.options()
#
## tests/sigtraceprogress.R
#
## Tests of *.sigtrace and *.progress
#
## $Revision: 1.4 $ $Date: 2018/11/02 00:53:45 $

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
  ## other code blocks
  e <- mad.progress(redwood, nsim=5)
  e <- mad.progress(redwood, nsim=19, alpha=0.05)
  f <- dclf.progress(redwood, nsim=5, scale=function(x) x^2)
  f <- dclf.progress(redwood, nsim=5, normalize=TRUE, deflate=TRUE)
  g <- dg.progress(redwood, nsim=5, scale=function(x) x^2)
  g <- dg.progress(redwood, nsim=5, normalize=TRUE, deflate=TRUE)
})
#'
#'     tests/simplepan.R
#'
#'   Tests of user interaction in simplepanel
#'   Handled by spatstatLocator()
#'
#'   $Revision: 1.2 $  $Date: 2018/10/16 00:46:41 $
#'

require(spatstat)

local({
  ## Adapted from example(simplepanel)
  ## make boxes
  outerbox <- owin(c(0,4), c(0,1))
  buttonboxes <- layout.boxes(outerbox, 4, horizontal=TRUE, aspect=1)
  ## make environment containing an integer count
  myenv <- new.env()
  assign("answer", 0, envir=myenv)
  ## what to do when finished: return the count.
  myexit <- function(e) { return(get("answer", envir=e)) }
  ## button clicks
  ## decrement the count
  Cminus <- function(e, xy) {
    ans <- get("answer", envir=e)
    assign("answer", ans - 1, envir=e)
    return(TRUE)
  }
  ## display the count (clicking does nothing)
  Cvalue <- function(...) { TRUE }
  ## increment the count
  Cplus <- function(e, xy) {
    ans <- get("answer", envir=e)
    assign("answer", ans + 1, envir=e)
    return(TRUE)
  }
  ## 'Clear' button
  Cclear <- function(e, xy) {
    assign("answer", 0, envir=e)
    return(TRUE)
  }
  ## quit button
  Cdone <- function(e, xy) { return(FALSE) }

  myclicks <- list("-"=Cminus,
                   value=Cvalue,
                   "+"=Cplus,
                   done=Cdone)
  ## redraw the button that displays the current value of the count
  Rvalue <- function(button, nam, e) {
     plot(button, add=TRUE)
     ans <- get("answer", envir=e)
     text(centroid.owin(button), labels=ans)
     return(TRUE)
  }
  ## make the panel
  P <- simplepanel("Counter",
                   B=outerbox, boxes=buttonboxes,
                   clicks=myclicks,
                   redraws = list(NULL, Rvalue, NULL, NULL),
                   exit=myexit, env=myenv)
  ## queue up a sequence of inputs
  boxcentres <- do.call(concatxy, unname(lapply(buttonboxes[c(3,3,1,3,2,4)],
                                                centroid.owin)))
  spatstat.utils::queueSpatstatLocator(boxcentres$x, boxcentres$y)
  ## go
  run.simplepanel(P)
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
#'  Basic tests of code in sparse3Darray.R and sparsecommon.R
#'  $Revision: 1.16 $ $Date: 2018/07/06 01:54:17 $

require(spatstat)
local({
  #' forming arrays

  #' creation by specifying nonzero elements
  M <- sparse3Darray(i=1:3, j=c(3,1,2), k=4:2,
                     x=runif(3), dims=rep(4, 3))
  #' duplicate entries
  Mn <- sparse3Darray(i=c(1,1,2), j=c(2,2,1), k=c(3,3,2),
                     x=runif(3), dims=rep(3, 3))
  #' cumulate entries in duplicate positions
  Ms <- sparse3Darray(i=c(1,1,2), j=c(2,2,1), k=c(3,3,2),
                      x=runif(3), dims=rep(3, 3), strict=TRUE)

  #' print method
  print(M)
  
  #' conversion of other data
  A <- array(c(1,3,0,0,0,0,0,4,0,2,0,5,
               0,0,1,0,0,0,1,0,0,0,1,0),
             dim=c(3,4,2))
  A1 <- A[,,1]
  A2 <- A[,,2]
  Z <- A[integer(0), , ]
  
  #' array to sparse array
  AA <- as.sparse3Darray(A) # positive extent
  ZZ <- as.sparse3Darray(Z) # zero extent
  #' list of matrices to sparse array
  AA <- as.sparse3Darray(list(A1, A2))
  #' matrix to sparse array
  AA1 <- as.sparse3Darray(A1)
  #' vector to sparse array
  A11 <- A[,1,1]
  AA11 <- as.sparse3Darray(A11)

  #' 
  dim(AA) <- dim(AA) + 1

  I1 <- SparseIndices(A1)
  I11 <- SparseIndices(A11)
  
  if(require(Matrix)) {
    #' sparse matrices from Matrix package
    A1 <- as(A1, "sparseMatrix")
    A2 <- as(A2, "sparseMatrix")
    A11 <- as(A11, "sparseVector")
    #' convert a list of sparse matrices to sparse array
    AA <- as.sparse3Darray(list(A1, A2))
    #' sparse matrix to sparse array
    AA1 <- as.sparse3Darray(A1)
    #' sparse vector to sparse array
    AA11 <- as.sparse3Darray(A11)

    #' internals 
    E1  <- SparseEntries(A1)
    I1  <- SparseIndices(A1)
    I11 <- SparseIndices(A11)
    df <- data.frame(i=c(1,3,5), j=3:1, k=rep(2, 3), x=runif(3))
    aa <- EntriesToSparse(df, NULL)
    bb <- EntriesToSparse(df, 7)
    cc <- EntriesToSparse(df, c(7, 4))
    dd <- EntriesToSparse(df, c(7, 4, 3))
  }

  BB <- evalSparse3Dentrywise(AA + AA/2)

  MM <- bind.sparse3Darray(M, M, along=1)
  MM <- bind.sparse3Darray(M, M, along=2)
})

    
local({

  if(require(Matrix)) {

    M <- sparse3Darray(i=1:4, j=sample(1:4, replace=TRUE),
                       k=c(1,2,1,2), x=1:4, dims=c(5,5,2))

    M

    dimnames(M) <- list(letters[1:5], LETTERS[1:5], c("yes", "no"))
    M
    
    U <- aperm(M, c(1,3,2))
    U

    #' tests of [.sparse3Darray
    M[ 3:4, , ]
    M[ 3:4, 2:4, ]
    M[ 4:3, 4:2, 1:2]
    M[, 3, ]
    M[, 3, , drop=FALSE]
    M[c(FALSE,TRUE,FALSE,FALSE,TRUE), , ]
    M[, , c(FALSE,FALSE), drop=FALSE]
    # matrix index
    M[cbind(3:5, 3:5, c(1,2,1))]
    M[cbind(3:5, 3:5, 2)]
    M[cbind(3:5,   2, 2)]
    M[cbind(c(2,2,4), c(3,3,2), 1)] # repeated indices
    
    MA <- as.array(M)
    UA <- as.array(U)

    ## tests of "[<-.sparse3Darray"
    Mflip <- Mzero <- MandM <- Mnew <- M
    Mflip[ , , 2:1] <- M
    stopifnot(Mflip[3,1,1] == M[3,1,2])
    Mzero[1:3,1:3,] <- 0
    stopifnot(all(Mzero[1,1,] == 0))
    M2a <- M[,,2,drop=FALSE]
    M2d <- M[,,2,drop=TRUE]
    MandM[,,1] <- M2a
    MandM[,,1] <- M2d
    ## slices of different dimensions
    M[ , 3, 1] <- 1:5
    M[2,  , 2] <- 1:5
    M[ 1, 3:5, 2] <- 4:6
    M[ 2, 5:3, 2] <- 4:6
    V3 <- sparseVector(x=1, i=2, length=3)
    M[ 1, 3:5, 2] <- V3
    M[ 2, 5:3, 2] <- V3
    M[,,2] <- M2a
    M[,,2] <- (M2a + 1)
    V5 <- sparseVector(x=1:2, i=2:3, length=5)
    M[,2,2] <- V5
    M[,,2] <- V5
    ## integer matrix index
    Mnew[cbind(3:5, 3:5, c(1,2,1))] <- 1:3
    Mnew[cbind(3:5, 3:5, 2)] <- 1:3
    Mnew[cbind(3:5,   2, 2)] <- 1:3
    Mnew[cbind(3:5, 3:5, c(1,2,1))] <- V3
    Mnew[cbind(3:5, 3:5, 2)] <- V3
    Mnew[cbind(3:5,   2, 2)] <- V3

    ## tests of arithmetic (Math, Ops, Summary)
    negM <- -M
    oneM <- 1 * M
    oneM <- M * 1
    twoM <- M + M
    range(M)

    cosM <- cos(M)  # non-sparse
    sinM <- sin(M)  # sparse

    Mpos <- (M > 0) # sparse
    Mzero <- !Mpos # non-sparse

    stopifnot(all((M+M) == 2*M))     # non-sparse
    stopifnot(!any((M+M) != 2*M))    # sparse

    ztimesM <- (1:5) * M  # sparse
    zplusM <- (1:5) + M  # non-sparse

    ## reconcile dimensions
    Msub <- M[,,1,drop=FALSE]
    Mdif <- M - Msub

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

    #' -----------  sparsecommon.R -----------------------
    B <- sparseMatrix(i=1:3, j=3:1, x= 10 * (1:3), dims=c(4,4))
    #' (and using sparse 3D array M and sparse vector V from above)

    Bmap <- mapSparseEntries(B, 1, 4:1)
    Mmap1 <- mapSparseEntries(M, 1, 5:1, across=3)
    Mmap2 <- mapSparseEntries(M, 3, 2:1, conform=FALSE)
    Mmap3 <- mapSparseEntries(M, 1, matrix(1:10, 5, 2), across=3)
    
    Vthrice  <- expandSparse(V, 3)
    VthriceT <- expandSparse(V, 3, 1)

    VV <- sparseVectorCumul(rep(1:3,2), rep(c(3,1,2), 2), 5)

    Vsum <- applySparseEntries(V, sum)
    Bdouble <- applySparseEntries(B, function(x) { 2 * x })
    Mminus <- applySparseEntries(M, function(x) -x)

    #'  -------------- sparselinalg.R -------------------------
    U <- aperm(M,c(3,1,2))  # 2 x 5 x 5
    w <- matrix(0, 5, 5)
    w[cbind(1:3,2:4)] <- 0.5
    w <- as(w, "sparseMatrix")
    UU <- sumsymouterSparse(U, w)
  }
})


#
#  tests/splitpea.R
#
#  Check behaviour of split.ppp etc
#
#  Thanks to Marcelino de la Cruz
#
#  $Revision: 1.12 $  $Date: 2019/01/18 01:58:29 $
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
Y <- split(X, "mineral")
print(Y)
print(summary(Y))
Y[c(TRUE,FALSE,TRUE)]
Y[1:2]
Y[3] <- Y[1]
})
#'
#'    tests/ssf.R
#'
#'   Tests of 'ssf' class
#'
#'   $Revision: 1.2 $ $Date: 2018/10/21 04:05:33 $
#'

require(spatstat)
local({
  Y <- cells[1:5]
  X <- rsyst(Window(Y), 5)
  Z <- runifpoint(3, Window(Y))
  f1 <- ssf(X, nncross(X,Y,what="dist"))
  f2 <- ssf(X, nncross(X,Y,what="dist", k=1:2))
  image(f1)
  g1 <- as.function(f1)
  g1(Z)
  g2 <- as.function(f2)
  g2(Z)
  plot(f1, style="contour")
  plot(f1, style="imagecontour")
  contour(f1)
  apply.ssf(f2, 1, sum)
  range(f1)
  min(f1)
  max(f1)
  integral(f1, weights=tile.areas(dirichlet(X)))
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
#'   tests/tessera.R
#'   Tessellation code, not elsewhere tested
#'   $Revision: 1.4 $ $Date: 2019/03/14 03:45:58 $
#'
require(spatstat)
local({
  W <- owin()
  Wsub <- square(0.5)
  X <- runifpoint(7, W)
  A <- dirichlet(X)
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
  flay(rotate, angle=pi/3)
  flay(affine, mat=matrix(c(1,2,0,1), 2, 2), vec=c(1,2))
  ## 
  unitname(B) <- c("metre", "metres")
  unitname(B)
  print(B)
  Bsub <- B[c(3,5,7)]
  print(Bsub)
  tilenames(H) <- letters[seq_along(tilenames(H))]
  #'
  Pe <- intersect.tess(A, Wsub)
  Pm <- intersect.tess(A, as.mask(Wsub))
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
#  $Revision: 1.4 $  $Date: 2019/01/02 08:27:51 $
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

#' other kinds of covariates
mod3 <- ppm(X ~ x + offset(y))
parres(mod3, "offset(y)")
Z <- distmap(runifpoint(3))
parres(mod3, Z)
mod4 <- ppm(X ~ sin(x), data=solist(B=Z))
parres(mod4, "sin(x)")
parres(mod4, "B")

})
#'
#'     tests/threedee.R
#'
#'     Tests of 3D code 
#'
#'      $Revision: 1.5 $ $Date: 2019/02/21 01:35:06 $
#'

require(spatstat)
local({
  X <- runifpoint3(30)
  Y <- runifpoint3(20)
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
  #' older code
  co <- coords(X)
  xx <- co$x
  yy <- co$y
  zz <- co$z
  gg1 <- g3engine(xx, yy, zz, correction="Hanisch G3")
  gg2 <- g3engine(xx, yy, zz, correction="minus sampling")
  ff1 <- f3engine(xx, yy, zz, correction="no")
  ff2 <- f3engine(xx, yy, zz, correction="minus sampling")
})
#'    tests/trigraph.R
#'
#'   Tests for C code in trigraf.c
#'   
#'  $Revision: 1.3 $  $Date: 2018/09/23 09:37:29 $
#'
require(spatstat)
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

reset.spatstat.options()


#
# tests/triplets.R
#
# test code for triplet interaction and associated summary function Tstat 
#
# $Revision: 1.6 $ $Date: 2018/07/02 15:51:26 $
#
require(spatstat)
local({
  fit <- ppm(redwood ~1, Triplets(0.1))
  fit
  suffstat(fit)
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
})
#
#  tests/undoc.R
#
#   $Revision: 1.6 $   $Date: 2018/10/19 09:45:52 $
#
#  Test undocumented hacks, etc

require(spatstat)
local({
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
#  $Revision: 1.10 $  $Date: 2019/02/10 07:12:06 $
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

  ## Model with zero-length coefficient vector
  lam <- intensity(X)
  f <- function(x,y) { rep(lam, length(x)) }
  model0 <- ppm(X ~ offset(log(f)) - 1)
  dd <- vcov(model0)
  cc <- vcov(model0, what="corr")

  ## Other weird stuff
  su <- suffloc(ppm(X ~ x))
})
#
# tests/windows.R
#
# Tests of owin geometry code
#
#  $Revision: 1.9 $  $Date: 2019/02/10 06:52:45 $

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
})

reset.spatstat.options()
##
## tests/xysegment.R
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
