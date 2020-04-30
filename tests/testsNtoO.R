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
# tests/NAinCov.R
#
# Testing the response to the presence of NA's in covariates
#
# $Revision: 1.7 $ $Date: 2020/04/30 05:23:52 $

if(FULLTEST) {
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
  #' quantile.ewcdf
  f <- ewcdf(runif(100), runif(100))
  qf <- quantile(f, probs=c(0.1, NA, 0.8))
  #' quantile.density
  f <- density(runif(100))
  qf <- quantile(f, probs=c(0.1, NA, 0.8))
})
}
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
# Also test nnorient()
#
#   $Revision: 1.33 $  $Date: 2020/04/30 05:23:52 $
#


local({
  eps <- sqrt(.Machine$double.eps)
  f <- function(mat,k) { apply(mat, 1, function(z,n) { sort(z)[n]  }, n=k+1) }
  g <- function(mat,k) { apply(mat, 1, function(z,n) { order(z)[n] }, n=k+1) }

  ## .......  Two dimensions ................
  if(ALWAYS) {
    X <- runifpoint(24)

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
  }

  if(FULLTEST) {
    a <- nndist(X, method="test")
    b <- nnwhich(X, method="test")
    a <- nndist(X, method="test", k=1:2)
    b <- nnwhich(X, method="test", k=1:2)
    a2 <- nndist(cells[1:3], k=1:3)
    b2 <- nnwhich(cells[1:3], k=1:3)
    a3 <- nndist(cells[1])
    b3 <- nnwhich(cells[1])
    m <- factor((1:npoints(X)) %% 2 == 0)
    a4 <- nndist.default(X, by=m, k=2)
    b4 <- nnwhich.default(X, by=m, k=2)
  }

  if(ALWAYS) {
    ## nncross.ppp without options
    Y <- runifpoint(30)
    Y <- Y[nndist(Y) > 0.02]
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

    ## nncross with sort on x
    nc <- nncross(X,Y, sortby="x")
    ncd <- nc$dist
    ncw <- nc$which
    if(any(abs(ncd - cdd) > eps))
      stop("nncross(sortby=x)$dist does not agree with apply(crossdist(), 1, min)")
    if(any(ncw != cdw))
      stop("nncross(sortby=x)$which does not agree with apply(crossdist(), 1, which.min)")

    ## nncross with data pre-sorted on x
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

    ## sanity check for nncross with k > 1
  ndw <- nncross(X, Y, k=1:4, what="which")
    if(any(is.na(ndw)))
      stop("NA's returned by nncross.ppp(k > 1, what='which')")
    nnc4 <- nncross(X, Y, k=1:4)
    iswhich <- (substr(colnames(nnc4), 1, nchar("which")) == "which")
    ndw <- nnc4[,iswhich]
    if(any(is.na(ndw)))
      stop("NA's returned by nncross.ppp(k > 1)$which")
  
    ## test of correctness for nncross with k > 1
    flipcells <- flipxy(cells)
    calcwhich <- nncross(cells, flipcells, k=1:4, what="which")
    truewhich <- t(apply(crossdist(cells,flipcells), 1, order))[,1:4]
    if(any(calcwhich != truewhich))
      stop("nncross(k > 1) gives wrong answer")
  }

  if(ALWAYS) {
    #' cover some C code blocks
    Z <- runifpoint(50)
    X <- Z[1:30]
    Y <- Z[20:50]
    iX <- 1:30
    iY <- 20:50
    Ndw <- nncross(X,Y, iX, iY, k=3)
    Nw  <- nncross(X,Y, iX, iY, k=3, what="which")
    Nd  <- nncross(X,Y, iX, iY, k=3, what="dist")

    ## special cases
    nndist(X[FALSE])
    nndist(X[1])
    nndist(X[1:3], k=4)
    nndist(X[1:3], k=1:4)
    nnwhich(X[FALSE])
    nnwhich(X[1])
    nnwhich(X[1:3], k=4)
    nnwhich(X[1:3], k=1:4)
    nncross(X[1:3], Y[FALSE])
    nncross(X[1:3], Y[1])
    nncross(X[1:3], Y[1:3], k=4)
    nncross(X[1:3], Y[1:3], k=1:4)
  }
  
  ## .......  Three dimensions ................

  if(ALWAYS) {
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

    ff <- function(mat,k) { apply(mat, 1, function(z,n) { sort(z)[n]  }, n=k) }
    gg <- function(mat,k) { apply(mat, 1, function(z,n) { order(z)[n] }, n=k) }

    Y <- runifpoint3(20)
    Y <- Y[nndist(Y) > 0.02]
    DXY <- crossdist(X,Y)
    a <- nncross(X,Y)
    a <- nncross(X,Y, what="dist")
    a <- nncross(X,Y, what="which")
    if(any(a != gg(DXY, 1)))
      stop("incorrect result from nncross.pp3(what='which')")
    a2 <- nncross(X,Y, k=2)
    a2 <- nncross(X,Y, what="dist", k=2)
    a2 <- nncross(X,Y, what="which", k=2)
    if(any(a2 != gg(DXY, 2)))
      stop("incorrect result from nncross.pp3(k=2, what='which')")
    iX <- 1:42
    iZ <- 30:42
    Z <- X[iZ]
    b <- nncross(X, Z, iX=iX, iY=iZ)
    b <- nncross(X, Z, iX=iX, iY=iZ, what="which")
    b <- nncross(X, Z, iX=iX, iY=iZ, what="dist")
    b2 <- nncross(X, Z, iX=iX, iY=iZ, k=2)
    b2 <- nncross(X, Z, iX=iX, iY=iZ, what="which", k=2)
    b2 <- nncross(X, Z, iX=iX, iY=iZ, what="dist", k=2)
    e1 <- nncross(X, Y[1:3], k=2:4)
    c1 <- nncross(X, Y, sortby="var")
    c2 <- nncross(X, Y, sortby="x")
    c3 <- nncross(X, Y, sortby="y")
    c4 <- nncross(X, Y, sortby="z")
    Xsort <- X[order(coords(X)$x)]
    c5 <- nncross(Xsort, Y, is.sorted.X=TRUE, sortby="x")
    Ysort <- Y[order(coords(Y)$x)]
    c6 <- nncross(Xsort, Ysort, is.sorted.X=TRUE, is.sorted.Y=TRUE, sortby="x")
  }

  if(FULLTEST) {
    ## special cases
    nndist(X[FALSE])
    nndist(X[1])
    nndist(X[1:3], k=4)
    nndist(X[1:3], k=1:4)
    nnwhich(X[FALSE])
    nnwhich(X[1])
    nnwhich(X[1:3], k=4)
    nnwhich(X[1:3], k=1:4)
    nncross(X[1:3], Y[FALSE])
    nncross(X[1:3], Y[1])
    nncross(X[1:3], Y[1:3], k=4)
    nncross(X[1:3], Y[1:3], k=1:4)
  }
  
  ## .......  m dimensions ................

  if(ALWAYS) {
    B <- boxx(c(0,1),c(0,1),c(0,1),c(0,1))
    X <- runifpointx(42, B)
    Y <- runifpointx(50, B)
    Y <- Y[nndist(Y) > 0.02]
    DXY <- crossdist(X,Y)
  
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

    a <- nncross(X,Y)
    ncd <- nncross(X,Y, what="dist")
    ncw <- nncross(X,Y, what="which")
    if(any(ncw != gg(DXY, 1)))
      stop("incorrect result from nncross.ppx(what='which')")
    a2 <- nncross(X,Y, k=2)
    ncd <- nncross(X,Y, what="dist", k=2)
    ncw <- nncross(X,Y, what="which", k=2)
    if(any(ncw != gg(DXY, 2)))
      stop("incorrect result from nncross.ppx(k=2, what='which')")
  }

  if(FULLTEST) {
    ## special cases
    nndist(X[FALSE])
    nndist(X[1])
    nndist(X[1:3], k=4)
    nndist(X[1:3], k=1:4)
    nnwhich(X[FALSE])
    nnwhich(X[1])
    nnwhich(X[1:3], k=4)
    nnwhich(X[1:3], k=1:4)
    nncross(X[1:3], Y[FALSE])
    nncross(X[1:3], Y[1])
    nncross(X[1:3], Y[1:3], k=4)
    nncross(X[1:3], Y[1:3], k=1:4)
  }

  if(ALWAYS) {
    ## test of agreement between nngrid.h and knngrid.h
    ##    dimyx=23 (found by trial-and-error) ensures that there are no ties 
    a <- as.matrix(nnmap(cells, what="which", dimyx=23))
    b <- as.matrix(nnmap(cells, what="which", dimyx=23, k=1:2)[[1]])
    if(any(a != b))
      stop("algorithms in nngrid.h and knngrid.h disagree")
    
    ## minnndist correctness
    X <- redwood3
    eps <- sqrt(.Machine$double.eps)
    mfast <- minnndist(X)
    mslow <- min(nndist(X))
    if(abs(mfast-mslow) > eps)
      stop("minnndist(X) disagrees with min(nndist(X))")

    ## maxnndist correctness
    mfast <- maxnndist(X)
    mslow <- max(nndist(X))
    if(abs(mfast-mslow) > eps)
      stop("maxnndist(X) disagrees with max(nndist(X))")
  }

  if(ALWAYS) {
    ## minnndist, maxnndist code blocks
    Y <- superimpose(amacrine, amacrine[10:20])
    a <- maxnndist(Y, positive=TRUE)
    u <- maxnndist(Y, positive=TRUE, by=marks(Y))
    b <- minnndist(Y, positive=TRUE)
    v <- minnndist(Y, positive=TRUE, by=marks(Y))

    ## nnmap code blocks
    A <- nnmap(cells[FALSE])
    A <- nnmap(cells, sortby="var")
    A <- nnmap(cells, sortby="x")
    A <- nnmap(cells, sortby="y")
    B <- nnmap(cells[1:3], k=4)
    B <- nnmap(cells[1:3], k=2:4)
    D <- nnmap(cells, outputarray=TRUE)
  }

  if(ALWAYS) {
    #' tests for has.close()
    #' (the default method uses nndist or pairdist, and can be trusted!)
    a <- has.close(redwood, 0.05)
    b <- has.close.default(redwood, 0.05)
    if(any(a != b)) stop("Incorrect result for has.close(X, r)")

    a <- has.close(redwood, 0.05, periodic=TRUE)
    a <- has.close.default(redwood, 0.05, periodic=TRUE)
    if(any(a != b)) stop("Incorrect result for has.close(X, r, periodic=TRUE)")

    Y <- split(amacrine)
    a <- with(Y, has.close(on, 0.05, off))
    b <- with(Y, has.close.default(on, 0.05, off))
    if(any(a != b)) stop("Incorrect result for has.close(X, r, Y)")

    a <- with(Y, has.close(on, 0.05, off, periodic=TRUE))
    b <- with(Y, has.close.default(on, 0.05, off, periodic=TRUE))
    if(any(a != b)) stop("Incorrect result for has.close(X, r, Y, periodic=TRUE)")
  }

  if(ALWAYS) {
    b <- bdist.pixels(letterR, style="coords")
    d <- bdist.pixels(letterR, dimyx=64, method="interpreted")
  }

  if(FULLTEST) {
    #' test nnorient
    nnorient(cells, domain=erosion(Window(cells), 0.1))
    #' degenerate case
    X <- cells[nndist(cells) > bdist.points(cells)]
    f <- nnorient(X)
    #' nnclean
    A <- nnclean(shapley, k=17, edge.correct=TRUE)
    B <- nnclean(runifpoint3(300), 3)
    #' stienen set
    #' bug when disc radius is zero
    Y <- unmark(humberside)[40:100] # contains duplicated points
    stienen(Y)
    Z <- stienenSet(Y)
    #' other cases
    U <- stienen(cells[1])
    V <- stienenSet(cells, edge=FALSE)

    ## nnfun.ppp
    f <- nnfun(cells)
    Z <- as.im(f)
    d <- domain(f)
    f <- nnfun(amacrine, value="mark")
    d <- domain(f)
    Z <- as.im(f)
    f <- nnfun(longleaf, value="mark")
    d <- domain(f)
    Z <- as.im(f)
  }

})

