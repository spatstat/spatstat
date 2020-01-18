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
##
##  tests/segments.R
##   Tests of psp class and related code
##                      [SEE ALSO: tests/xysegment.R]
##
##  $Revision: 1.22 $  $Date: 2020/01/18 02:29:37 $

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

#' cases of superimpose.psp
A <- as.psp(matrix(runif(40), 10, 4), window=owin())
B <- as.psp(matrix(runif(40), 10, 4), window=owin())
superimpose(A, B, W=ripras)
superimpose(A, B, W="convex")

#' tests of density.psp
Y <- as.psp(simplenet)
YC <- density(Y, 0.2, method="C", edge=FALSE, dimyx=64)
YI <- density(Y, 0.2, method="interpreted", edge=FALSE, dimyx=64)
YF <- density(Y, 0.2, method="FFT", edge=FALSE, dimyx=64)
xCI <- max(abs(YC/YI - 1))
xFI <- max(abs(YF/YI - 1))
if(xCI > 0.01) stop(paste("density.psp C algorithm relative error =", xCI))
if(xFI > 0.01) stop(paste("density.psp FFT algorithm relative error =", xFI))
B <- square(0.3)
density(Y, 0.2, at=B)
density(Y, 0.2, at=B, edge=TRUE, method="C")
Z <- runifpoint(3, B)
density(Y, 0.2, at=Z)
density(Y, 0.2, at=Z, edge=TRUE, method="C")

#' as.psp.data.frame

  df <- as.data.frame(matrix(runif(40), ncol=4))
  A <- as.psp(df, window=square(1))
  colnames(df) <- c("x0","y0","x1","y1")
  df <- cbind(df, data.frame(marks=1:nrow(df)))
  B <- as.psp(df, window=square(1))
  colnames(df) <- c("xmid", "ymid", "length", "angle", "marks")
  E <- as.psp(df, window=square(c(-1,2)))
  G <- E %mark% factor(sample(letters[1:3], nsegments(E), replace=TRUE))
  H <- E %mark% runif(nsegments(E))
  
#' print and summary methods
  A
  B
  E
  G
  H
  summary(B)
  summary(G)
  summary(H)
  M <- B
  marks(M) <- data.frame(id=marks(B), len=lengths.psp(B))
  M
  summary(M)
  subset(M, select=len)

  #' plot method cases  
  spatstat.options(monochrome=TRUE)
  plot(B)
  plot(G)
  plot(M)
  spatstat.options(monochrome=FALSE)
  plot(B)
  plot(G)
  plot(M)
  #' misuse of 'col' argument - several cases
  plot(G, col="grey") # discrete
  plot(B, col="grey") 
  plot(unmark(B), col="grey") 
  plot(M, col="grey") 
  
  #' miscellaneous class support cases
  marks(M) <- marks(M)[1,,drop=FALSE]

  #' undocumented  
  as.ppp(B)

  #' segment crossing code
  X <- psp(runif(30),runif(30),runif(30),runif(30), window=owin())
  A <- selfcut.psp(X, eps=1e-11)
  B <- selfcut.psp(X[1])
  #' 
  Y <- psp(runif(30),runif(30),runif(30),runif(30), window=owin())
  Z <- edges(letterR)[c(FALSE,TRUE)]
  spatstat.options(selfcrossing.psp.useCall=FALSE, crossing.psp.useCall=FALSE)
  A <- selfcrossing.psp(X)
  B <- selfcrossing.psp(Z)
  D <- crossing.psp(X,Y,details=TRUE)
  spatstat.options(selfcrossing.psp.useCall=TRUE, crossing.psp.useCall=TRUE)
  A <- selfcrossing.psp(X)
  B <- selfcrossing.psp(Z)
  D <- crossing.psp(X,Y,details=TRUE)
})

local({
  #' test rshift.psp and append.psp with marks (Ute Hahn)
  m <- data.frame(A=1:10, B=letters[1:10])
  g <- gl(3, 3, length=10)
  X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin(), marks=m)
  Y <- rshift(X, radius = 0.1)
  Y <- rshift(X, radius = 0.1, group=g)
  #' mark management
  b <- data.frame(A=1:10)
  X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin(), marks=b)
  stopifnot(is.data.frame(marks(X)))
  Y <- rshift(X, radius = 0.1)
  Y <- rshift(X, radius = 0.1, group=g)
})

local({
  #' geometry
  m <- data.frame(A=1:10, B=letters[1:10])
  X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin(), marks=m)
  Z <- rotate(X, angle=pi/3, centre=c(0.5, 0.5))
  Y <- endpoints.psp(X, which="lower")
  Y <- endpoints.psp(X, which="upper")
  Y <- endpoints.psp(X, which="right")
  U <- flipxy(X)
})

local({
  ## nnfun.psp
  P <- psp(runif(10), runif(10), runif(10), runif(10),
           window=square(1), marks=runif(10))
  f <- nnfun(P)
  f <- nnfun(P, value="mark")
  d <- domain(f)
  Z <- as.im(f)
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
# $Revision: 1.2 $ $Date: 2020/01/10 04:54:49 $
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
  extractAIC(fit)
  fitx <- update(fit, . ~ x)
  simulate(fitx, seed=42)
  unitname(fitx)
  unitname(fitx) <- "km"

  mur <- solapply(murchison,rescale, 1000, "km")
  mur$dfault <- distfun(mur$faults)
  fut <- slrm(gold ~ dfault, data=mur, splitby="greenstone")
  A <- model.images(fut)
})


#'    tests/sparse3Darrays.R
#'  Basic tests of code in sparse3Darray.R and sparsecommon.R
#'  $Revision: 1.21 $ $Date: 2019/12/31 02:38:48 $

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
  #' NULL with warning
  as.sparse3Darray(list())

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
    #' duplicated entries
    dfdup <- df[c(1:3, 2), ]
    aa <- EntriesToSparse(dfdup, NULL)
    bb <- EntriesToSparse(dfdup, 7)
    cc <- EntriesToSparse(dfdup, c(7, 4))
    dd <- EntriesToSparse(dfdup, c(7, 4, 3))
  }

  BB <- evalSparse3Dentrywise(AA + AA/2)

  MM <- bind.sparse3Darray(M, M, along=1)
  MM <- bind.sparse3Darray(M, M, along=2)

  RelevantEmpty(42)
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
    M[1:2, 1, 2:3] # exceeds array bounds
    # matrix index
    M[cbind(3:5, 3:5, c(1,2,1))]
    M[cbind(3:5, 3:5, 2)]
    M[cbind(3:5,   2, 2)]
    M[cbind(c(2,2,4), c(3,3,2), 1)] # repeated indices
    M[cbind(1:4, 1, 2:3)] # exceeds array bounds

    MA <- as.array(M)
    UA <- as.array(U)

    Mfix <- sparse3Darray(i=1:3, j=c(3,1,2), k=4:2,
                          x=runif(3), dims=rep(4, 3))
    Mfix[cbind(1,3,4)] # single entry - occupied
    Mfix[cbind(1,2,4)] # single entry - unoccupied
    Mfix[cbind(1,c(2,3,2,3),4)] # sparse vector with repeated entries
    

    ## tests of "[<-.sparse3Darray"
    Mflip <- Mzero <- MandM <- Mnew <- Mext <- M
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
    Mext[1,2,3] <- 4 # exceeds array bounds
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
    Mduf <- Msub - M
    
    ## tensor operator
    tenseur(c(1,-1), M, 1, 3)
    tenseur(M, M, 1:2, 1:2)
    tenseur(M, M, 1:2, 2:1)
    V <- sparseVector(i=c(1,3,6),x=1:3, length=7)
    tenseur(V,V)
    tenseur(V,V,1,1)
    A <- sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=c(7, 15))
    A[1:4, 2:5] <- 3
    tenseur(A, A, 1, 1)
    tenseur(t(A), A, 2, 1)
    tenseur(V, A, 1, 1)
    tenseur(t(A), V, 2, 1)
    tenseur(as.vector(V), A, 1, 1)
    tenseur(t(A), as.vector(V), 2, 1)

    v <- 0:3
    tensor1x1(v, Mfix)
    tensor1x1(as(v, "sparseVector"), Mfix)
    
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

    ## complex valued arrays
    Mcplx <- sparse3Darray(i=1:3, j=c(3,1,2), k=4:2,
                           x=runif(3)+runif(3)*1i, dims=rep(4, 3))
    print(Mcplx)
    

    #' -----------  sparsecommon.R -----------------------
    B <- sparseMatrix(i=1:3, j=3:1, x= 10 * (1:3), dims=c(4,4))
    #' (and using sparse 3D array M and sparse vector V from above)
    V2 <- sparseVector(i=c(2,3,6),x=4:6, length=7)  # different pattern
    check.anySparseVector(V2, 10, fatal=FALSE)

    Bmap <- mapSparseEntries(B, 1, 4:1)
    Mmap1 <- mapSparseEntries(M, 1, 5:1, across=3)
    Mmap2 <- mapSparseEntries(M, 3, 2:1, conform=FALSE)
    Mmap3 <- mapSparseEntries(M, 1, matrix(1:10, 5, 2), across=3)
    
    Vmap <- mapSparseEntries(V, 1, V2)
    Vmap <- mapSparseEntries(V, 1, 8)
    Vthrice  <- expandSparse(V, 3)
    VthriceT <- expandSparse(V, 3, 1)
    VF <- as.vector(V) # non-sparse
    VFmap <- mapSparseEntries(VF, 1, V2)
    VFmap <- mapSparseEntries(VF, 1, 8)
    VFthrice  <- expandSparse(VF, 3)
    VFthriceT <- expandSparse(VF, 3, 1)
    VFthriceX <- expandSparse(VF, 3, 2)
    
    VV <- sparseVectorCumul(rep(1:3,2), rep(c(3,1,2), 2), 5)

    Vsum <- applySparseEntries(V, sum)
    Bdouble <- applySparseEntries(B, function(x) { 2 * x })
    Mminus <- applySparseEntries(M, function(x) -x)

    # empty sparse matrices/arrays
    Bempty <- B
    Bempty[] <- 0
    mapSparseEntries(Bempty, 1, 42)
    Mempty <- M
    Mempty[] <- 0
    Mmap1 <- mapSparseEntries(Mempty, 1, 5:1, across=3)
    Mmap2 <- mapSparseEntries(Mempty, 3, 2:1, conform=FALSE)
    Mmap3 <- mapSparseEntries(Mempty, 1, matrix(1:10, 5, 2), across=3)

    #'  -------------- sparselinalg.R -------------------------
    U <- aperm(M,c(3,1,2))  # 2 x 5 x 5
    w <- matrix(0, 5, 5)
    w[cbind(1:3,2:4)] <- 0.5
    w <- as(w, "sparseMatrix")
    UU <- sumsymouterSparse(U, w, dbg=TRUE)
    Uempty <- sparse3Darray(dims=c(2,5,5))
    UU <- sumsymouterSparse(Uempty, w, dbg=TRUE)
  }
})

local({
  # 1 x 1 x 1 arrays
  M1 <- sparse3Darray(i=1, j=1, k=1, x=42, dims=rep(1,3))
  M0 <- sparse3Darray(                     dims=rep(1,3))
  i1 <- matrix(1, 1, 3)
  a1 <- M1[i1]
  a0 <- M0[i1]
  A <- array(runif(75) * (runif(75) < 0.7), dim=c(3,5,5))
  M <- as.sparse3Darray(A)
  M[rep(1,3), c(1,1,2), rep(2, 3)]
})
#
#  tests/splitpea.R
#
#  Check behaviour of split.ppp etc
#
#  Thanks to Marcelino de la Cruz
#
#  $Revision: 1.13 $  $Date: 2019/12/15 04:46:57 $
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

# other bugs/ code blocks in split.ppp, split<-.ppp, [<-.splitppp
flog <- rep(c(TRUE,FALSE), 21)
fimg <- as.im(dirichlet(runifpoint(5, Window(cells))))
A <- split(cells, flog)
B <- split(cells, square(0.5))
D <- split(cells, fimg)
E <- split(cells, logical(42), drop=TRUE)
Cellules <- cells
split(Cellules, flog) <- solapply(A, rjitter)
split(Cellules, fimg) <- solapply(D, rjitter)
D[[2]] <- rjitter(D[[2]])
Funpines <- finpines
marks(Funpines)[,"diameter"] <- factor(marks(Funpines)[,"diameter"])
G <- split(Funpines)
H <- split(Funpines, "diameter")
split(Funpines) <- solapply(G, rjitter)
split(Funpines, "diameter") <- solapply(H, rjitter)

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


#'
#'   tests/sumfun.R
#'
#'   Tests of code for summary functions
#'       including score residual functions etc
#'
#'   $Revision: 1.3 $ $Date: 2020/01/05 01:55:10 $

require(spatstat)
local({
  W <- owin(c(0,1), c(-1/2, 0))
  Gr <- Gest(redwood, correction="all",domain=W)
  Fr <- Fest(redwood, correction="all",domain=W)
  Jr <- Jest(redwood, correction="all",domain=W)
  
  F0 <- Fest(redwood[FALSE], correction="all")
  Fh <- Fest(humberside, domain=erosion(Window(humberside), 100))

  FIr <- Finhom(redwood, savelambda=TRUE)
  JIr <- Jinhom(redwood, savelambda=TRUE)
  
  Ga <- Gcross(amacrine, correction="all")
  Ia <- Iest(amacrine, correction="all")
  lam <- intensity(amacrine)
  lmin <- 0.9 * min(lam)
  nJ <- sum(marks(amacrine) == "off")
  FM <- FmultiInhom(amacrine, marks(amacrine) == "off",
                    lambdaJ=rep(lam["off"], nJ),
                    lambdamin = lmin)
  GM <- GmultiInhom(amacrine, marks(amacrine) == "on",
                    marks(amacrine) == "off",
                    lambda=lam[marks(amacrine)],
                    lambdamin=lmin,
                    ReferenceMeasureMarkSetI=42)

  pt <- psst(cells, interaction=Strauss(0.1), fun=nndcumfun)

  a <- compileCDF(D=nndist(redwood),
                  B=bdist.points(redwood),
                  r=seq(0, 1, length=256))
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
  X <- runifpoint(8)

  ## symbolmap 
  g1 <- symbolmap(range=c(0,100), size=function(x) x/50)
  invoke.symbolmap(g1, 50, x=numeric(0), y=numeric(0), add=TRUE)
  plot(g1, labelmap=100)
  ## constant/trivial
  a <- symbolmap(pch=16)
  print(a)
  plot(a)
  symbolmapdomain(a)
  b <- symbolmap()
  print(b)

  ## textureplot
  V <- as.im(dirichlet(X))
  tmap <- textureplot(V)
  textureplot(V, textures=tmap, legend=TRUE, leg.side="left")
  textureplot(V, leg.side="bottom")
  textureplot(V, leg.side="top")
  ## spacing too large for tiles - upsets various pieces of code
  textureplot(V, spacing=2)

  ## plot.texturemap
  plot(tmap, vertical=TRUE)
  plot(tmap, vertical=TRUE, xlim=c(0,1))
  plot(tmap, vertical=TRUE, ylim=c(0,1))
  plot(tmap, vertical=FALSE, xlim=c(0,1))
  plot(tmap, vertical=FALSE, ylim=c(0,1))

  ## infrastructure
  plan.legend.layout(owin(), side="top", started=TRUE)
})
