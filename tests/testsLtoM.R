#'
#'   tests/layered.R
#'
#'   Tests of 'layered' class
#'
#'   $Revision: 1.1 $  $Date: 2018/07/14 06:23:45 $
#'
require(spatstat)
local({
  D <- distmap(cells)
  L <- layered(D, cells,
               plotargs=list(list(ribbon=FALSE), list(pch=16)))
  #'
  plot(L, which=2, plotargs=list(list(pch=3)))
  plot(L, plotargs=list(list(pch=3)))
  #'
  W <- as.owin(L)
  V <- domain(L)
  #' methods
  L2 <- L[square(0.5)]
  Lr <- reflect(L)
  Lf <- flipxy(L)
  Ls <- scalardilate(L, 2)
  La <- shift(L, origin="midpoint")
  Lo <- rotate(L, pi/3, origin="bottomleft")
  Lu <- rescale(L, 0.1, "parsec")
  #' as.layered 
  M <- as.layered(finpines)
  M2 <- as.layered(split(amacrine))
})
## 
##    tests/legacy.R
##
## Test that current version of spatstat is compatible with outmoded usage
## $Revision: 1.2 $ $Date: 2015/12/29 08:54:49 $
local({
  require(spatstat)

  ## (1) Old syntax of ppm
  ppm(cells, ~x)
  
  ## (2) Old syntax of MultiStrauss etc.
  r <- matrix(3, 2, 2)
  a <- MultiStrauss( , r)
  a <- MultiStrauss(NULL, r)
  a <- MultiHard(, r)
  
  h <- r/2
  a <- MultiStraussHard( , r, h)

  NULL
})
#'
#'    tests/leverinf.R
#'
#'   leverage and influence for Gibbs models
#' 
#'   $Revision: 1.24 $ $Date: 2019/02/11 09:40:50 $
#' 

require(spatstat)
local({
  cat("Running non-sparse algorithm...", fill=TRUE)
  # original non-sparse algorithm
  Leverage <- function(...) leverage(..., sparseOK=FALSE)
  Influence <- function(...) influence(..., sparseOK=FALSE)
  Dfbetas <- function(...) dfbetas(..., sparseOK=FALSE)
  # Strauss()$delta2
  fitS <- ppm(cells ~ x, Strauss(0.12), rbord=0)
  levS <- Leverage(fitS)
  infS <- Influence(fitS)
  dfbS <- Dfbetas(fitS)
  # Geyer()$delta2
  fitG <- ppm(redwood ~ 1, Geyer(0.1, 2), rbord=0)
  levG <- Leverage(fitG)
  infG <- Influence(fitG)
  # AreaInter()$delta2
  fitA <- ppm(cells ~ 1, AreaInter(0.07), rbord=0, nd=11)
  levA <- Leverage(fitA)
  infA <- Influence(fitA)
  # pairwise.family$delta2
  fitD <- ppm(cells ~ 1, DiggleGatesStibbard(0.12), rbord=0)
  levD <- Leverage(fitD)
  infD <- Influence(fitD)
  # DiggleGratton() special code
  fitDG <- ppm(cells ~ 1, DiggleGratton(0.05, 0.12), rbord=0)
  levDG <- Leverage(fitDG)
  infDG <- Influence(fitDG)
  # ppmInfluence; offset is present; coefficient vector has length 0
  fitH <- ppm(cells ~ 1, Hardcore(0.07))
  levH <- Leverage(fitH)
  infH <- Influence(fitH)
  # ppmInfluence; hard core
  fitSH <- ppm(cells ~ 1, StraussHard(0.07, 0.01))
  levSH <- Leverage(fitSH)
  infSH <- Influence(fitSH)
  # ppmInfluence; offset is present; coefficient vector has length 1
  fitHx <- ppm(cells ~ x, Hardcore(0.07), rbord=0)
  levHx <- Leverage(fitHx)
  infHx <- Influence(fitHx)
  ## multitype 
  futAm <- ppm(amacrine ~ x + marks, Strauss(0.07))
  levAm <- leverage(futAm)

  ## .........   class support .............................
  ## other methods for classes leverage.ppm and influence.ppm
  ## not elsewhere tested
  cat("Testing class support...", fill=TRUE)
  w <- domain(levS)
  w <- Window(infS)
  vv <- shift(levS, c(1.2, 1.3))
  vv <- shift(infS, c(1.2, 1.3))
  A <- quadrats(Window(cells), 2)
  a <- integral(levS,domain=A)
  b <- integral(infS,domain=A)
  u <- Smooth(levS, sigma=0.07)
  v <- Smooth(infS, sigma=0.1)
  ## plot options
  plot(levS, what="exact")
  plot(levS, what="nearest")
  contour(levS, what="nearest")
  persp(levS, what="nearest")
  ## plotting for multitype models
  plot(levAm)
  contour(levAm)
  persp(levAm)
  plot(levAm, multiplot=FALSE)
  contour(levAm, multiplot=FALSE)

  ## ..........  compare algorithms .........................
  ## divide and recombine algorithm
  cat("Reduce maximum block side to 50,000 ...", fill=TRUE)
  op <- spatstat.options(maxmatrix=50000)
  ## non-sparse
  levSB <- Leverage(fitS)
  infSB <- Influence(fitS)
  dfbSB <- Dfbetas(fitS)

  chk <- function(x, y, what,
                  from="single-block and multi-block",
                  thresh=1e-12) {
    if(max(abs(x-y)) > thresh)
      stop(paste("Different results for", what, "obtained from",
                 from, "algorithms"),
           call.=FALSE)
    invisible(NULL)
  }

  cat("Compare single-block to multi-block...", fill=TRUE)
  chk(marks(as.ppp(infS)), marks(as.ppp(infSB)), "influence")
  chk(as.im(levS),         as.im(levSB),         "leverage")
  chk(dfbS$val,            dfbSB$val,            "dfbetas$value")
  chk(dfbS$density,        dfbSB$density,        "dfbetas$density")

  # also check case of zero cif
  cat("Check zero cif cases...", fill=TRUE)
  levHB <- Leverage(fitH)
  infHB <- Influence(fitH)
  dfbHB <- Dfbetas(fitH)
  levHxB <- Leverage(fitHx)
  infHxB <- Influence(fitHx)
  dfbHxB <- Dfbetas(fitHx)

  ## run all code segments
  Everything <- function(model, ...) { ppmInfluence(model, ..., what="all") }

  cat("Run full code on AreaInteraction model...", fill=TRUE)
  pmiA <- Everything(fitA)
  
  ## sparse algorithm, with blocks
  cat("Run sparse algorithm with blocks...", fill=TRUE)
  pmiSSB <- Everything(fitS, sparseOK=TRUE)
  # also check case of zero cif
  pmiHSB <- Everything(fitH, sparseOK=TRUE)
  pmiSHSB <- Everything(fitSH, sparseOK=TRUE)
  pmiHxSB <- Everything(fitHx, sparseOK=TRUE)

  cat("Reinstate maxmatrix...", fill=TRUE)
  spatstat.options(op)

  ## sparse algorithm, no blocks
  cat("Compare sparse and non-sparse results...", fill=TRUE)
  pmi <- Everything(fitS, sparseOK=TRUE)
  levSp <- pmi$leverage
  infSp <- pmi$influence
  dfbSp <- pmi$dfbetas
  chks <- function(...) chk(..., from="sparse and non-sparse")
  
  chks(marks(as.ppp(infS)), marks(as.ppp(infSp)), "influence")
  chks(as.im(levS),         as.im(levSp),         "leverage")
  chks(dfbS$val,            dfbSp$val,            "dfbetas$value")
  chks(dfbS$density,        dfbSp$density,        "dfbetas$density")

  #' case of zero cif
  cat("zero cif...", fill=TRUE)
  pmiH <- Everything(fitH, sparseOK=TRUE)
  pmiSH <- Everything(fitSH, sparseOK=TRUE)
  pmiHx <- Everything(fitHx, sparseOK=TRUE)

  #' other code blocks - check execution only
  cat("other code blocks...", fill=TRUE)
  a <- Everything(fitS) 
  a <- Everything(fitS, method="interpreted") 
  a <- Everything(fitS, method="interpreted", entrywise=FALSE)
  a <- Everything(fitS,                       entrywise=FALSE)
  #' zero cif
  b <- Everything(fitSH) 
  b <- Everything(fitSH, method="interpreted") 
  b <- Everything(fitSH, method="interpreted", entrywise=FALSE)
  b <- Everything(fitSH,                       entrywise=FALSE) 
  #' NOTE: code for irregular parameters is tested below, and in 'make bookcheck'

  ## ...........  logistic fits .......................
  cat("Logistic fits...", fill=TRUE)
  #'  special algorithm for delta2
  fitSlogi <- ppm(cells ~ x, Strauss(0.12), rbord=0, method="logi")
  pmiSlogi <- Everything(fitSlogi)
  #'  special algorithm for delta2
  fitGlogi <- ppm(redwood ~ 1, Geyer(0.1, 2), rbord=0, method="logi")
  pmiGlogi <- Everything(fitGlogi)
  #'  generic algorithm for delta2
  fitDlogi <- ppm(cells ~ 1, DiggleGatesStibbard(0.12), rbord=0, method="logi")
  pmiDlogi <- Everything(fitDlogi)
  #'  generic algorithm for delta2 : offset; zero-dimensional 
  fitHlogi <- ppm(cells ~ 1, Hardcore(0.07), method="logi")
  pmiHlogi <- Everything(fitHlogi)
  #'  generic algorithm for delta2 : offset; 1-dimensional 
  fitHxlogi <- ppm(cells ~ x, Hardcore(0.07), rbord=0, method="logi")
  pmiHxlogi <- Everything(fitHxlogi)
  #' plotting
  plot(leverage(fitSlogi))
  plot(influence(fitSlogi))
  plot(dfbetas(fitSlogi))
  
  #' other code blocks - check execution only
  cat("Other code blocks...", fill=TRUE)
  b <- Everything(fitSlogi)  # i.e. full set of results
  b <- Everything(fitSlogi, method="interpreted") 
  b <- Everything(fitSlogi, method="interpreted", entrywise=FALSE)
  b <- Everything(fitSlogi,                       entrywise=FALSE) 

  #' irregular parameters
  cat("Irregular parameters...", fill=TRUE)
  ytoa <- function(x,y, alpha=1) { y^alpha }
  lam <- function(x,y,alpha=1) { exp(4 + y^alpha) }
  set.seed(90210)
  X <- rpoispp(lam, alpha=2)
  iScor <- list(alpha=function(x,y,alpha) { alpha * y^(alpha-1) } )
  iHess <- list(alpha=function(x,y,alpha) { alpha * (alpha-1) * y^(alpha-2) } )
  gogo <- function(tag, ..., iS=iScor, iH=iHess) {
    #' compute full set of results
    cat(tag, fill=TRUE)
    ppmInfluence(..., what="all", iScore=iS, iHessian=iH)
  }
  cat("Offset model...", fill=TRUE)
  fut <- ippm(X ~ offset(ytoa), start=list(alpha=1))
  d <- gogo("a", fut)
  d <- gogo("b", fut, method="interpreted") 
  d <- gogo("c", fut, method="interpreted", entrywise=FALSE)
  d <- gogo("d", fut,                       entrywise=FALSE) 
  cat("Offset+x model...", fill=TRUE)
  futx <- ippm(X ~ x + offset(ytoa), start=list(alpha=1))
  d <- gogo("a", futx) 
  d <- gogo("b", futx, method="interpreted") 
  d <- gogo("c", futx, method="interpreted", entrywise=FALSE)
  d <- gogo("d", futx,                       entrywise=FALSE)
  cat("Offset model Strauss ...", fill=TRUE)
  futS <- ippm(X ~ offset(ytoa), Strauss(0.07), start=list(alpha=1))
  d <- gogo("a", futS)
  d <- gogo("b", futS, method="interpreted") 
  d <- gogo("c", futS, method="interpreted", entrywise=FALSE)
  d <- gogo("d", futS,                       entrywise=FALSE) 
  cat("Offset+x model Strauss ...", fill=TRUE)
  futxS <- ippm(X ~ x + offset(ytoa), Strauss(0.07), start=list(alpha=1))
  d <- gogo("a", futxS) 
  d <- gogo("b", futxS, method="interpreted") 
  d <- gogo("c", futxS, method="interpreted", entrywise=FALSE)
  d <- gogo("d", futxS,                       entrywise=FALSE)

  #'
  set.seed(452)
  foo <- ppm(cells ~ 1, Strauss(0.15), method="ho", nsim=5)
  aa <- Everything(foo)
})

reset.spatstat.options()
##
##    tests/linalgeb.R
##
## checks validity of linear algebra code
##
##  $Revision: 1.4 $ $Date: 2019/02/21 02:21:43 $
##
require(spatstat)
local({
  p <- 3
  n <- 4
  k <- 2
  
  x <- matrix(1:(n*p), n, p)
  w <- runif(n)
  y <- matrix(1:(2*n), n, k)
  zUS <- zWS <- matrix(0, p, p)
  zUA <- zWA <- matrix(0, p, k)
  for(i in 1:n) {
    zUS <- zUS +        outer(x[i,],x[i,])
    zWS <- zWS + w[i] * outer(x[i,],x[i,])
    zUA <- zUA +        outer(x[i,],y[i,])
    zWA <- zWA + w[i] * outer(x[i,],y[i,])
  }
  if(!identical(zUS, sumouter(x)))
    stop("sumouter gives incorrect result in Unweighted Symmetric case")
  if(!identical(zWS, sumouter(x,w)))
    stop("sumouter gives incorrect result in Weighted Symmetric case")
  if(!identical(zUA, sumouter(x, y=y)))
    stop("sumouter gives incorrect result in Unweighted Asymmetric case")
  if(!identical(zWA, sumouter(x, w, y)))
    stop("sumouter gives incorrect result in Weighted Asymmetric case")
  
  x <- array(as.numeric(1:(p * n * n)), dim=c(p, n, n))
  w <- matrix(1:(n*n), n, n)
  y <- matrix(numeric(p * p), p, p)
  for(i in 1:n)
    for(j in (1:n)[-i])
      y <- y + w[i,j] * outer(x[,i,j], x[,j,i])
  z <- sumsymouter(x, w)
  if(!identical(y,z))
    stop("sumsymouter gives incorrect result")
})
## 
## tests/localpcf.R
##
## temporary test file for localpcfmatrix
##  $Revision: 1.2 $  $Date: 2015/12/29 08:54:49 $

require(spatstat)
local({
  a <- localpcfmatrix(redwood)
  a
  plot(a)
  a[, 3:5]
})
#
# tests/lppstuff.R
#
# Tests for lpp code
#
#  $Revision: 1.38 $  $Date: 2019/10/29 04:14:33 $


require(spatstat)

local({
  #' lpp class support
  X <- runiflpp(5, simplenet)
  Xone <- X %mark% runif(5)
  Xtwo <- X %mark% data.frame(a=1:5, b=runif(5))
  print(summary(Xone))
  print(summary(Xtwo))
  plot(X, show.window=TRUE)
  plot(Xone)
  plot(Xtwo, do.several=FALSE)
  #' geometry etc
  rotate(X, pi/3, centre=c(0.2,0.3))
  superimpose.lpp(L=simplenet)
  #' cut.lpp
  tes <- lineardirichlet(runiflpp(4, simplenet))
  f <- as.linfun(tes)
  Z <- as.linim(f)
  cut(X, tes)
  cut(X, f)
  cut(X, Z)

  #' check 'normalise' option in linearKinhom
  fit <- lppm(X ~x)
  K <- linearKinhom(X, lambda=fit, normalise=FALSE)
  plot(K)
  g <- linearpcfinhom(X, lambda=fit, normalise=FALSE)
  plot(g)
  K <- linearKinhom(X, lambda=fit, normalise=TRUE)
  plot(K)
  g <- linearpcfinhom(X, lambda=fit, normalise=TRUE)
  plot(g)
  ## other code blocks
  K <- linearKinhom(X, lambda=fit, correction="none", ratio=TRUE)
  K0 <- linearKcross(dendrite[1], "thin", "thin")
  # check empty patterns OK
  X <- runiflpp(0, simplenet)
  print(X)
  
  ## nearest neighbour distances
  eps <- sqrt(.Machine$double.eps)
  f <- function(mat,k) { apply(mat, 1, function(z,n) { sort(z)[n]  }, n=k+1) }
  g <- function(mat,k) { apply(mat, 1, function(z,n) { order(z)[n] }, n=k+1) }

  XX <- spiders
  nn <- nndist(XX)
  nnP <- f(pairdist(XX), 1)
  if(any(abs(nn - nnP) > eps))
    stop("nndist.lpp does not agree with pairdist.lpp")

  nw <- nnwhich(XX)
  nwP <- g(pairdist(XX), 1)
  if(any(nw != nwP))
    stop("nnwhich.lpp does not agree with pairdist")

  #' code blocks in nndist.lpp/nnwhich.lpp
  #' non-sparse network, interpreted code  
  Ad <- nndist(spiders, method="interpreted") 
  Aw <- nnwhich(spiders, method="interpreted")
  #' sparse network, older C code
  opa <- spatstat.options(Cnndistlpp=FALSE)
  Bd <- nndist(dendrite) 
  Bw <- nnwhich(dendrite) 
  spatstat.options(opa)
  #' undefined nearest neighbours
  Ed <- nndist(spiders[1:3], k=1:3)
  Ew <- nnwhich(spiders[1:3], k=1:3)
    
  #' trivial cases in nncross.lpp
  a <- nncross(runiflpp(0, simplenet), runiflpp(1, simplenet),
               what="which", format="list")$which
  a <- nncross(runiflpp(0, simplenet), runiflpp(1, simplenet),
               what="dist", format="list")$dist

  #' compare algorithms             
  ZZ <- split(chicago)
  XX <- ZZ$damage
  YY <- ZZ$assault
  op <- spatstat.options(Cnncrosslpp=FALSE)
  a <- nncross(XX, YY)
  spatstat.options(Cnncrosslpp=TRUE)
  b <- nncross(XX, YY)
  if(any(a$which != b$which))
    stop("Inconsistent values of nncross.lpp()$which from different C code")
  if(max(abs(a$dist - b$dist)) > eps)
    stop("Inconsistent values of nncross.lpp()$dist from different C code")

  spatstat.options(Cnncrosslpp=TRUE)
  b2 <- nncross(XX, YY, k=1:2, what="which")
  if(any(b2$which.1 != b$which))
    stop("inconsistent values of nncross.lpp()$which from k=1:2 and k=1")
  a2 <- nncross(XX, YY, k=1:2, what="dist")
  if(max(abs(a2$dist.1 - a$dist)) > eps)
    stop("Inconsistent values of nncross.lpp()$dist from k=1:2 and k=1")

  spatstat.options(Cnncrosslpp=TRUE)
  ii <- seq_len(npoints(XX))
  w1 <- nnwhich(XX)
  w2 <- nncross(XX, XX, iX=ii, iY=ii, what="which")
  w3 <- nncross(XX, XX, iX=ii, iY=ii, what="which", method="interpreted")
  if(any(w1 != w2))
    stop("nnwhich.lpp disagrees with nncross.lpp(iX, iY)")
  if(any(w2 != w3))
    stop("Different results for nncross.lpp(iX, iY, 'which') using R and C")
  d1 <- nndist(XX)
  d2 <- nncross(XX, XX, iX=ii, iY=ii, what="dist")
  d3 <- nncross(XX, XX, iX=ii, iY=ii, what="dist", method="interpreted")
  if(max(abs(d1-d2)) > eps)
    stop("nndist.lpp disagrees with nncross.lpp(iX, iY)")
  if(max(abs(d2-d3)) > eps)
    stop("Different results for nncross.lpp(iX, iY, 'dist') using R and C")
  spatstat.options(Cnncrosslpp=FALSE)
  w4 <- nncross(XX, XX, iX=ii, iY=ii, what="which")
  d4 <- nncross(XX, XX, iX=ii, iY=ii, what="dist")
  if(any(w2 != w4))
    stop("Different results for nncross.lpp(iX, iY, 'which') fast and slow C")
  if(max(abs(d2-d4)) > eps)
    stop("Different results for nncross.lpp(iX, iY, 'dist') fast and slow C")
  spatstat.options(Cnncrosslpp=TRUE)

  spatstat.options(op)
  reset.spatstat.options()

  # test handling marginal cases
  xyd <- nncross(XX, YY[1])
  A <- runiflpp(5, simplenet)
  B <- runiflpp(2, simplenet)
  aaa <- nncross(A,B,k=3:5) #' all undefined
  spatstat.options(Cnncrosslpp=FALSE)
  bbb <- nncross(B,A, iX=1:2, iY=1:5) # another code block
  spatstat.options(Cnncrosslpp=TRUE)

  reset.spatstat.options()
  
  ## as.linnet.psp (Suman's example)
  Lines <- as.data.frame(as.psp(simplenet))
  newseg <- c(Lines[1,1:2], Lines[10,3:4])
  Lines <- rbind(Lines, newseg)
  Y <- as.psp(Lines, window=Window(simplenet))
  marks(Y) <- c(3, 4, 5, 5, 3, 4, 5, 5,5, 5,1)
  Z <- as.linnet(Y) # can crash if marks don't match segments
  
  ## Test linnet surgery code
  set.seed(42)
  X <- runiflpp(30, simplenet)
  V <- runiflpp(30, simplenet)
  XV <- insertVertices(X, V)
  validate.lpp.coords(XV, context="calculated by insertVertices")

  ## Test [.lpp internal data
  B <- owin(c(0.1,0.7),c(0.19,0.6))
  XB <- X[B]
  validate.lpp.coords(XB, context="returned by [.lpp")

  ## Tests related to linearK, etc
  testcountends <- function(X, r=100, s=1) {
    if(s != 1) {
      X <- rescale(X, s)
      r <- r/s
    }
    L <- as.linnet(X)
    n1 <- countends(L, X[1], r)
    n2 <- npoints(lineardisc(L, X[1], r, plotit=FALSE)$endpoints)
    if(n1 != n2)
      stop(paste("Incorrect result from countends:",
                 n1, "!=", n2, 
                 paren(paste("scale=", 1/s))),
           call.=FALSE)
  }
  # original scale
  X <- unmark(chicago)
  testcountends(X)
  # finer scale
  testcountends(X, s=1000)
  #' disconnected
  L <- thinNetwork(simplenet, retainedges = -c(3,8))
  S <- as.psp(L)
  x <- midpoints.psp(S)[1]
  len <- lengths.psp(S)[1]
  A <- lineardisc(L, x, len,  plotit=FALSE) # involves many segments of network
  B <- lineardisc(L, x, len/5, plotit=FALSE) # involves one segment of network
  op <- spatstat.options(Ccountends=FALSE)
  A <- lineardisc(L, x, len,  plotit=FALSE)
  B <- lineardisc(L, x, len/5, plotit=FALSE)
  spatstat.options(op)
  reset.spatstat.options()
  
  ## Test algorithms for boundingradius.linnet
  L <- as.linnet(chicago, sparse=TRUE)
  L$boundingradius <- NULL # artificially remove
  opa <- spatstat.options(Clinearradius=FALSE)
  bR <- as.linnet(L, sparse=FALSE)$boundingradius
  spatstat.options(Clinearradius=TRUE)
  bC <- as.linnet(L, sparse=FALSE)$boundingradius
  spatstat.options(opa)
  if(abs(bR-bC) > 0.001 * (bR+bC)/2)
    stop("Disagreement between R and C algorithms for boundingradius.linnet",
         call.=FALSE)

  ## linnet things
  is.connected(as.linnet(dendrite))
  zik <- rescale(chicago, 39.37/12, "m")
  Simon <- simplenet
  unitname(Simon) <- list("metre", "metres", 0.5)
  b <- rescale(Simon)
  ds <- density(simplenet, 0.05)

  ## integral.linim with missing entries
  xcoord <- linfun(function(x,y,seg,tp) { x }, domain(chicago))
  xcoord <- as.linim(xcoord, dimyx=32)
  integral(xcoord)

  ## options to plot.linim
  plot(xcoord, legend=FALSE)
  plot(xcoord, leg.side="top")
  
  ## as.linim.linim
  xxcc <- as.linim(xcoord)
  xxcceps <- as.linim(xcoord, eps=15)
  xxccdel <- as.linim(xcoord, delta=30)
  df1 <- attr(xxcc, "df")
  df2 <- attr(xxccdel, "df")
  df3 <- resampleNetworkDataFrame(df1, df2)

  ## linim with complex values
  Zc <- as.im(function(x,y){(x-y) + x * 1i}, Frame(simplenet))
  Fc <- linim(simplenet, Zc)
  print(Fc)
  print(summary(Fc))
  
  ## linim with df provided
  Z <- as.im(function(x,y) {x-y}, Frame(simplenet))
  X <- linim(simplenet, Z)
  df <- attr(X, "df")
  XX <- linim(simplenet, Z, df=df)
  dfwithout <- df[, colnames(df) != "values"]
  XXX <- linim(simplenet, Z, df=dfwithout)
  plot(XXX, zlim=c(-1,1))
  plot(XXX, legend=FALSE)
  plot(XXX, leg.side="bottom")
  
  ## lpp with multiple columns of marks
  M <- chicago
  marks(M) <- cbind(type=marks(M), data.frame(distnearest=nndist(M)))
  plot(M, main="")
  summary(M)
  MM <- cut(M)

  #' other cases
  CX <- cut(chicago)
  nd <- nndist(spiders)
  SX <- cut(spiders %mark% nd, breaks=3)
  SX <- cut(spiders, nd, breaks=c(0,100,200,Inf), include.lowest=TRUE)
  
  ## linequad
  X <- runiflpp(6, simplenet)
  Y <- X %mark% factor(rep(c("A", "B"), 3))
  aX <- linequad(X)
  aY <- linequad(Y)
  aXR <- linequad(X, random=TRUE)
  aYR <- linequad(Y, random=TRUE)
  P <- as.ppp(X)
  S <- as.psp(domain(X))
  d <- linequad(P, S)
  oop <- spatstat.options(Clinequad=FALSE)
  bX <- linequad(X)
  spatstat.options(oop)

  ## other internal utilities
  df <- pointsAlongNetwork(simplenet, 0.05)
  X <- as.ppp(df[,c("x", "y")], W=Frame(simplenet))
  A <- local2lpp(simplenet, seg=df$seg, tp=df$tp, X=X, df.only=FALSE)

  ## mark-mark scatterplot uses pairdist
  X <- runiflpp(20, simplenet) %mark% runif(20)
  markmarkscatter(X, 0.2)
  markmarkscatter(X[FALSE], 0.1)

  ## tree branches
  # make a simple tree
  m <- simplenet$m
  m[8,10] <- m[10,8] <- FALSE
  L <- linnet(vertices(simplenet), m)
  tb <- treebranchlabels(L, 1)
  X <- runiflpp(50, L)
  # delete branch B
  XminusB <- deletebranch(X, "b", tb)
  # extract branch B
  XB <- extractbranch(X, "b", tb)

  ## cases of lintess()
  A <- lintess(simplenet) # argument df missing 
  S <- as.psp(simplenet)
  ns <- nsegments(S)
  df <- data.frame(seg=1:ns, t0=0, t1=1, tile=letters[1:ns])
  M <- data.frame(len=lengths.psp(S), ang=angles.psp(S))
  V <- lintess(simplenet, df, marks=M)

  ## methods for class lintess
  U <- unmark(V)
  U <- unstack(V)
  print(summary(V))
  W <- Window(V)
  plot(V, style="image")
  plot(V, style="width")

  ## linear tessellations infrastructure
  nX <- 100
  nY <- 20
  X <- runiflpp(nX, simplenet)
  Y <- runiflpp(nY, simplenet)
  tes <- divide.linnet(Y)
  cX <- coords(X)
  iI <- lineartileindex(cX$seg, cX$tp, tes, method="interpreted")
  iC <- lineartileindex(cX$seg, cX$tp, tes, method="C")
  iE <- lineartileindex(cX$seg, cX$tp, tes, method="encode")
  if(!identical(iI,iC))
    stop("conflicting results from lineartileindex (interpreted vs C)")
  if(!identical(iI,iE))
    stop("conflicting results from lineartileindex (interpreted vs encoded)")
  iA <- as.linfun(tes)(X)
  if(!identical(iI, iA))
    stop("disagreement between as.linfun.lintess and lineartileindex")
})

reset.spatstat.options()

local({
  #' handling by 'solist', 'unstack', 'plot.solist' etc
  L <- simplenet
  X <- runiflpp(5, L) %mark% cbind(a=1:5, b=letters[1:5])
  ns <- nsegments(L)
  df <- data.frame(seg=1:ns, t0=0, t1=1, tile=letters[1:ns])
  S <- lintess(L, df)
  f <- as.linfun(S)
  g <- as.linfun(S, values=seq_len(nsegments(L)))
  V <- as.linim(f)
  Z <- as.linim(g)
  shebang <- solist(L=L, X=X, S=S, f=f, g=g, V=V, Z=Z)
  plot(shebang)
  plot(shebang, valuesAreColours=FALSE)
  kapow <- unstack(shebang)
  plot(kapow)
})

reset.spatstat.options()

local({
  #' densityVoronoi.lpp and related code
  X <- runiflpp(5, simplenet)
  densityVoronoi(X, f=0)
  densityVoronoi(X, f=1e-8)
  densityVoronoi(X, f=1)
  densityVoronoi(X[FALSE], f=0.5)
  XX <- X[rep(1:5, 4)]
  densityVoronoi(XX, f=0.99999, nrep=5)
  #' bandwidth selection
  bw.voronoi(X, nrep=4, prob=c(0.2, 0.4, 0.6))
})

local({
  #' complex-valued functions and images
  f <- function(x,y,seg,tp) { x + y * 1i }
  g <- linfun(f, simplenet)
  h <- as.linim(g)
  plot(Re(h))
  plot(h)
  plot(g)
  integral(h)
  integral(g)
})
  
local({
  ## bug in 'lixellate' (Jakob Gulddahl Rasmussen)
  X <- ppp(c(0,1), c(0,0), owin())
  L <- linnet(X, edges = matrix(1:2, ncol=2))
  Y <- lpp(X, L)
  ## The left end point is OK
  lixellate(Y[1], nsplit=30)
  d <- density(Y[1], .1)
  ## The right end point gave an error
  lixellate(Y[2], nsplit=30)
  d <- density(Y[2], .1)
})

local({
  ## make some bad data and repair it
  X <- runiflpp(4, simplenet)
  sx1 <- coords(X)$seg[1]
  ns <- nsegments(X)
  X$domain$from <- X$domain$from[c(1:ns, sx1)]
  X$domain$to   <- X$domain$to[c(1:ns, sx1)]
  Y <- repairNetwork(X)
})
#'
#'   lppmodels.R
#'
#'   Tests of lppm and class support
#' 
#'   $Revision: 1.1 $ $Date: 2018/05/13 04:14:28 $
#'

require(spatstat)

local({
  fit0 <- lppm(spiders)
  fit1 <- lppm(spiders ~ x)
  fit2 <- lppm(chicago ~ x+y)
  X <- runiflpp(10, simplenet)
  Z <- distfun(runiflpp(10, simplenet))
  fit3 <- lppm(X ~ Z)

  summary(fit0)
  summary(fit1)
  summary(fit2)
  summary(fit3)
  
  pseudoR2(fit0)
  pseudoR2(fit1)
  pseudoR2(fit2)
  pseudoR2(fit3)

  Window(fit1)

  a <- model.images(fit0)
  a <- model.images(fit1)
  a <- model.images(fit2)
  a <- model.images(fit3)

  b <- model.matrix(fit0)
  b <- model.matrix(fit1)
  b <- model.matrix(fit2)
  b <- model.matrix(fit3)

  is.multitype(fit0)
  is.multitype(fit1)
  is.multitype(fit2)
  is.multitype(fit3)

  fit0e <- emend(fit0)
  fit1e <- emend(fit1)
  fit2e <- emend(fit2)
  fit3e <- emend(fit3)

  #' fundamental utilities:
  #' evalCovar
  ycoord <- function(x,y) { y }
  YS <- as.linim(ycoord, L=domain(spiders))
  YC <- as.linim(ycoord, L=domain(chicago))

  aT <- evalCovar(fit1, YS, interpolate=TRUE)
  aF <- evalCovar(fit1, YS, interpolate=FALSE)
  dT <- evalCovar(fit1, ycoord, interpolate=TRUE)
  dF <- evalCovar(fit1, ycoord, interpolate=FALSE)

  bT <- evalCovar(fit2, YC, interpolate=TRUE)
  bF <- evalCovar(fit2, YC, interpolate=FALSE)
  cT <- evalCovar(fit2, ycoord, interpolate=TRUE)
  cF <- evalCovar(fit2, ycoord, interpolate=FALSE)
})
##
##     tests/marcelino.R
##
##     $Revision: 1.3 $  $Date: 2015/12/29 08:54:49 $
##
require(spatstat)

local({
  Y <- split(urkiola)
  B <- Y$birch
  O <- Y$oak
  B.lam <- predict (ppm(B ~polynom(x,y,2)), type="trend")
  O.lam <- predict (ppm(O ~polynom(x,y,2)), type="trend")

  Kinhom(B, lambda=B.lam, correction="iso")
  Kinhom(B, lambda=B.lam, correction="border")

  Kcross.inhom(urkiola, i="birch", j="oak", B.lam, O.lam)
  Kcross.inhom(urkiola, i="birch", j="oak", B.lam, O.lam, correction = "iso")
  Kcross.inhom(urkiola, i="birch", j="oak", B.lam, O.lam, correction = "border")
})


##
##    tests/markcor.R
##
##   Tests of mark correlation code (etc)
##
## $Revision: 1.5 $ $Date: 2019/08/01 06:32:20 $

require(spatstat)

local({
  ## check.testfun checks equality of functions
  ##  and is liable to break if the behaviour of all.equal is changed
  fe <- function(m1, m2) {m1 == m2}
  fm <- function(m1, m2) {m1 * m2}
  fs <- function(m1, m2) {sqrt(m1)}
  if(check.testfun(fe, X=amacrine)$ftype != "equ")
    warning("check.testfun fails to recognise mark equality function")
  if(check.testfun(fm, X=longleaf)$ftype != "mul")
    warning("check.testfun fails to recognise mark product function")
  check.testfun(fs, X=longleaf)
  check.testfun("mul")
  check.testfun("equ")
  
  ## test all is well in Kmark -> Kinhom 
  MA <- Kmark(amacrine,function(m1,m2){m1==m2})
  set.seed(42)
  AR <- rlabel(amacrine)
  MR <- Kmark(AR,function(m1,m2){m1==m2})
  if(isTRUE(all.equal(MA,MR)))
    stop("Kmark unexpectedly ignores marks")

  ## cover code blocks in markcorr()
  X <- runifpoint(100) %mark% runif(100)
  Y <- X %mark% data.frame(u=runif(100), v=runif(100))
  ww <- runif(100)
  fone <- function(x) { x/2 }
  ffff <- function(x,y) { fone(x) * fone(y) }
  aa <- markcorr(Y)
  bb <- markcorr(Y, ffff, weights=ww, normalise=TRUE)
  bb <- markcorr(Y, ffff, weights=ww, normalise=FALSE)
  bb <- markcorr(Y, f1=fone, weights=ww, normalise=TRUE)
  bb <- markcorr(Y, f1=fone, weights=ww, normalise=FALSE)

  ## markcrosscorr
  a <- markcrosscorr(betacells, normalise=FALSE)
  if(require(sm)) {
    b <- markcrosscorr(betacells, method="sm")
  }
})
#' tests/mctests.R
#' Monte Carlo tests
#'        (mad.test, dclf.test, envelopeTest, hasenvelope)
#' $Revision: 1.2 $ $Date: 2019/06/03 10:39:31 $

require(spatstat)
local({
  envelopeTest(cells, Lest, exponent=1, nsim=9, savepatterns=TRUE)
  (a3 <- envelopeTest(cells, Lest, exponent=3, nsim=9, savepatterns=TRUE))

  envelopeTest(a3, Lest, exponent=3, nsim=9, alternative="less")
  
  fitx <- ppm(redwood~x)
  envelopeTest(fitx, exponent=2, nsim=9, savefuns=TRUE)

  envelopeTest(redwood, Lest, exponent=1, nsim=19,
               rinterval=c(0, 0.1), alternative="greater", clamp=TRUE)
  envelopeTest(redwood, pcf, exponent=Inf, nsim=19,
               rinterval=c(0, 0.1), alternative="greater", clamp=TRUE)
})

#'  tests/morpho.R
#' 
#' morphology code blocks
#'
#' $Revision: 1.2 $ $Date: 2018/05/13 04:02:40 $

require(spatstat)
local({
  #' owin
  a <- erosion(letterR, 0.1, polygonal=FALSE)
  b <- dilation(letterR, 0.1, polygonal=FALSE)
  at <- erosion(letterR, 0.1, polygonal=FALSE, strict=TRUE)
  bt <- dilation(letterR, 0.1, polygonal=FALSE, tight=FALSE)
  #' psp
  S <- edges(letterR)
  dm <- dilation(S, 0.1, polygonal=FALSE)
  dt <- dilation(S, 0.1, polygonal=FALSE, tight=FALSE)
  op <- spatstat.options(old.morpho.psp=TRUE)
  dn <- dilation(S, 0.1, polygonal=TRUE)
  spatstat.options(op)
  cS <- closing(S, 0.1, polygonal=FALSE)
  eS <- erosion(S, 0)
  oS <- opening(S, 0)
  #' ppp
  dc <- dilation(cells, 0.06, polygonal=FALSE)
  ec <- erosion(cells, 0)
  oc <- opening(cells, 0)
})

reset.spatstat.options()
#
# tests/mppm.R
#
# Basic tests of mppm
#
# $Revision: 1.13 $ $Date: 2019/12/06 10:01:08 $
# 

require(spatstat)

local({
  ## test interaction formulae and subfits
  fit1 <- mppm(Points ~ group, simba,
               hyperframe(po=Poisson(), str=Strauss(0.1)),
               iformula=~ifelse(group=="control", po, str))
  fit2 <- mppm(Points ~ group, simba,
               hyperframe(po=Poisson(), str=Strauss(0.1)),
               iformula=~str/id)
  fit2w <- mppm(Points ~ group, simba,
                hyperframe(po=Poisson(), str=Strauss(0.1)),
                iformula=~str/id, weights=runif(nrow(simba)))
# currently invalid  
#  fit3 <- mppm(Points ~ group, simba,
#               hyperframe(po=Poisson(), pie=PairPiece(c(0.05,0.1))),
#        iformula=~I((group=="control") * po) + I((group=="treatment") * pie))
  fit1
  fit2
  fit2w
#  fit3

  ## run summary.mppm which currently sits in spatstat-internal.Rd
  summary(fit1)
  summary(fit2)
  summary(fit2w)
#  summary(fit3)

  ## test vcov algorithm
  vcov(fit1)
  vcov(fit2)
#  vcov(fit3)

  ## test subfits algorithm
  s1 <- subfits(fit1)
  s2 <- subfits(fit2)
#  s3 <- subfits(fit3)

  ## validity of results of subfits()
  p1 <- solapply(s1, predict)
  p2 <- solapply(s2, predict)
#  p3 <- solapply(s3, predict)

})

local({
  ## cases of predict.mppm
  W <- solapply(waterstriders, Window)
  Fakes <- solapply(W, runifpoint, n=30)
  FakeDist <- solapply(Fakes, distfun)
  H <- hyperframe(Bugs=waterstriders,
                  D=FakeDist)
  fit <- mppm(Bugs ~ D, data=H)
  p1 <- predict(fit)
  p2 <- predict(fit, locations=Fakes)
  p3 <- predict(fit, locations=solapply(W, erosion, r=4))
  locn <- as.data.frame(do.call(cbind, lapply(Fakes, coords)))
  df <- data.frame(id=sample(1:3, nrow(locn), replace=TRUE),
                   D=runif(nrow(locn)))
  p4 <- predict(fit, locations=locn, newdata=df)
})

local({
  ##  [thanks to Sven Wagner]
  ## factor covariate, with some levels unused in some rows
  set.seed(14921788)
  H <- hyperframe(X=replicate(3, runifpoint(20), simplify=FALSE),
                  Z=solist(as.im(function(x,y){x}, owin()),
                    as.im(function(x,y){y}, owin()),
                    as.im(function(x,y){x+y}, owin())))
  H$Z <- solapply(H$Z, cut, breaks=(0:4)/2)

  fit6 <- mppm(X ~ Z, H)
  v6 <- vcov(fit6)
  s6 <- subfits(fit6)
  p6 <- solapply(s6, predict)

  # random effects
  fit7 <- mppm(X ~ Z, H, random=~1|id)
  v7 <- vcov(fit7)
  s7 <- subfits(fit7)
  p7 <- solapply(s7, predict)

  fit7a <- mppm(X ~ Z, H, random=~x|id)
  v7a <- vcov(fit7a)
  s7a <- subfits(fit7a)
  p7a <- solapply(s7a, predict)

  # multitype: collisions in vcov.ppm, predict.ppm
  H$X <- lapply(H$X, rlabel, labels=factor(c("a","b")), permute=FALSE)
  M <- MultiStrauss(matrix(0.1, 2, 2), c("a","b"))
  fit8 <- mppm(X ~ Z, H, M)
  v8 <- vcov(fit8, fine=TRUE)
  s8 <- subfits(fit8)
  p8 <- lapply(s8, predict)
  c8 <- lapply(s8, predict, type="cif")

  fit9 <- mppm(X ~ Z, H, M, iformula=~Interaction * id)
  v9 <- vcov(fit9, fine=TRUE)
  s9 <- subfits(fit9)
  p9 <- lapply(s9, predict)
  c9 <- lapply(s9, predict, type="cif")

  # and a simple error in recognising 'marks'
  fit10 <- mppm(X ~ marks, H)
})

local({
  ## test handling of offsets and zero cif values in mppm
  H <- hyperframe(Y = waterstriders)
  mppm(Y ~ 1,  data=H, Hardcore(1.5))
  mppm(Y ~ 1,  data=H, StraussHard(7, 1.5))

  ## prediction, in training/testing context
  ##    (example from Markus Herrmann and Ege Rubak)
  X <- waterstriders
  dist <- solapply(waterstriders,
                   function(z) distfun(runifpoint(1, Window(z))))
  i <- 3
  train <- hyperframe(pattern = X[-i], dist = dist[-i])
  test <- hyperframe(pattern = X[i], dist = dist[i])
  fit <- mppm(pattern ~ dist, data = train)
  pred <- predict(fit, type="cif", newdata=test, verbose=TRUE)
})

local({
  ## test handling of interaction coefficients in multitype case
  set.seed(42)
  XX <- as.solist(replicate(3, rthin(amacrine, 0.8), simplify=FALSE))
  H <- hyperframe(X=XX)
  M <- MultiStrauss(matrix(0.1, 2, 2), levels(marks(amacrine)))
  fit <- mppm(X ~ 1, H, M)
  co <- coef(fit)
  subco <- sapply(subfits(fit), coef)
  if(max(abs(subco - co)) > 0.001)
    stop("Wrong coefficient values in subfits, for multitype interaction")
})

local({
  ## test lurking.mppm
  # example from 'mppm'
  n <- 7
  H <- hyperframe(V=1:n,
                  U=runif(n, min=-1, max=1))
  H$Z <- setcov(square(1))
  H$U <- with(H, as.im(U, as.rectangle(Z)))
  H$Y <- with(H, rpoispp(eval.im(exp(2+3*Z))))
  fit <- mppm(Y ~ Z + U + V, data=H)

  lurking(fit, expression(Z), type="P")
  lurking(fit, expression(V), type="raw") # design covariate
  lurking(fit, expression(U), type="raw") # image, constant in each row
  lurking(fit, H$Z,           type="P")   # list of images
})

local({
  ## test anova.mppm code blocks and scoping problem
  H <- hyperframe(X=waterstriders)
  mod0 <- mppm(X~1, data=H, Poisson())
  modxy <- mppm(X~x+y, data=H, Poisson())
  mod0S <- mppm(X~1, data=H, Strauss(2))
  modxyS <- mppm(X~x+y, data=H, Strauss(2))
  anova(mod0, modxy, test="Chi")
  anova(mod0S, modxyS, test="Chi")
  anova(modxy, test="Chi")
  anova(modxyS, test="Chi")
})

local({
  ## test multitype stuff
  foo <- flu[1:3,]
  msh <- MultiStraussHard(iradii=matrix(100, 2, 2),
                          hradii=matrix(10,2,2),
                          types=levels(marks(foo$pattern[[1]])))
  msh0 <- MultiStraussHard(iradii=matrix(100, 2, 2),
                          hradii=matrix(10,2,2))
  fit <- mppm(pattern ~ 1, data=foo, interaction=msh0)
  print(fit)
  print(summary(fit))
})
#'
#'     tests/msr.R
#'
#'     $Revision: 1.1 $ $Date: 2019/02/15 01:37:55 $
#'
#'     Tests of code for measures
#'

require(spatstat)

local({
  rr <- residuals(ppm(cells ~ x))

  a <- summary(rr)
  b <- is.marked(rr)
  w <- as.owin(rr)
  z <- domain(rr)
  ss <- scalardilate(rr, 2)

  rrr <- augment.msr(rr, sigma=0.08)
  uuu <- update(rrr)

  mm <- residuals(ppm(amacrine ~ x))
  ss <- residuals(ppm(amacrine ~ x), type="score")

  plot(mm)
  plot(mm, multiplot=FALSE)
  plot(mm, equal.markscale=TRUE, equal.ribbon=TRUE)
  plot(ss)
  plot(ss, multiplot=FALSE)
})
