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
#'   tests/layered.R
#'
#'   Tests of 'layered' class
#'
#'   $Revision: 1.2 $  $Date: 2020/04/29 08:55:17 $
#'
if(FULLTEST) {
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
}
## 
##    tests/legacy.R
##
## Test that current version of spatstat is compatible with outmoded usage
## $Revision: 1.3 $ $Date: 2020/04/29 08:55:17 $

if(FULLTEST) {
local({

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
}
#'
#'    tests/leverinf.R
#'
#'   leverage and influence for Gibbs models
#' 
#'   $Revision: 1.29 $ $Date: 2021/01/22 08:09:13 $
#' 

if(!FULLTEST)
  spatstat.options(npixel=32, ndummy.min=16)

local({
  cat("Running non-sparse algorithm...", fill=TRUE)
  # original non-sparse algorithm
  Leverage <- function(...) leverage(..., sparseOK=FALSE)
  Influence <- function(...) influence(..., sparseOK=FALSE)
  Dfbetas <- function(...) dfbetas(..., sparseOK=FALSE)
  if(ALWAYS) {
    ## Strauss()$delta2
    fitS <- ppm(cells ~ x, Strauss(0.12), rbord=0)
    levS <- Leverage(fitS)
    infS <- Influence(fitS)
    dfbS <- Dfbetas(fitS)
    ## Geyer()$delta2
    fitG <- ppm(redwood ~ 1, Geyer(0.1, 2), rbord=0)
    levG <- Leverage(fitG)
    infG <- Influence(fitG)
    ## AreaInter()$delta2
    fitA <- ppm(cells ~ 1, AreaInter(0.06), rbord=0, nd=11)
    levA <- Leverage(fitA)
    infA <- Influence(fitA)
    ## pairwise.family$delta2
    fitD <- ppm(cells ~ 1, DiggleGatesStibbard(0.12), rbord=0)
    levD <- Leverage(fitD)
    infD <- Influence(fitD)
    ## DiggleGratton() special code
    fitDG <- ppm(cells ~ 1, DiggleGratton(0.05, 0.12), rbord=0)
    levDG <- Leverage(fitDG)
    infDG <- Influence(fitDG)
    ## ppmInfluence; offset is present; coefficient vector has length 0
    fitH <- ppm(cells ~ 1, Hardcore(0.07))
    levH <- Leverage(fitH)
    infH <- Influence(fitH)
    ## ppmInfluence; hard core
    fitSH <- ppm(cells ~ 1, StraussHard(0.07, 0.01))
    levSH <- Leverage(fitSH)
    infSH <- Influence(fitSH)
    ## ppmInfluence; offset is present; coefficient vector has length 1
    fitHx <- ppm(cells ~ x, Hardcore(0.07), rbord=0)
    levHx <- Leverage(fitHx)
    infHx <- Influence(fitHx)
    ## multitype 
    futAm <- ppm(amacrine ~ x + marks, Strauss(0.07))
    levAm <- leverage(futAm)
  }

  if(FULLTEST) {
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
  }

  if(ALWAYS) {
    ## ..........  compare algorithms .........................
    ## divide and recombine algorithm
    cat("Reduce maximum block side to 50,000 ...", fill=TRUE)
    op <- spatstat.options(maxmatrix=50000)
    ## non-sparse
    levSB <- Leverage(fitS)
    infSB <- Influence(fitS)
    dfbSB <- Dfbetas(fitS)
  }

  chk <- function(x, y, what,
                  from="single-block and multi-block",
                  thresh=1e-12) {
    if(max(abs(x-y)) > thresh)
      stop(paste("Different results for", what, "obtained from",
                 from, "algorithms"),
           call.=FALSE)
    invisible(NULL)
  }

  if(ALWAYS) {
    cat("Compare single-block to multi-block...", fill=TRUE)
    chk(marks(as.ppp(infS)), marks(as.ppp(infSB)), "influence")
    chk(as.im(levS),         as.im(levSB),         "leverage")
    chk(dfbS$val,            dfbSB$val,            "dfbetas$value")
    chk(dfbS$density,        dfbSB$density,        "dfbetas$density")

    ## also check case of zero cif
    cat("Check zero cif cases...", fill=TRUE)
    levHB <- Leverage(fitH)
    infHB <- Influence(fitH)
    dfbHB <- Dfbetas(fitH)
    levHxB <- Leverage(fitHx)
    infHxB <- Influence(fitHx)
    dfbHxB <- Dfbetas(fitHx)
  }

  ## run all code segments
  Everything <- function(model, ...) { ppmInfluence(model, ..., what="all") }
    
  if(FULLTEST) {
    cat("Run full code on AreaInteraction model...", fill=TRUE)
    pmiA <- Everything(fitA)
    
    ## sparse algorithm, with blocks
    cat("Run sparse algorithm with blocks...", fill=TRUE)
    pmiSSB <- Everything(fitS, sparseOK=TRUE)
    ## also check case of zero cif
    pmiHSB <- Everything(fitH, sparseOK=TRUE)
    pmiSHSB <- Everything(fitSH, sparseOK=TRUE)
    pmiHxSB <- Everything(fitHx, sparseOK=TRUE)

    cat("Reinstate maxmatrix...", fill=TRUE)
    spatstat.options(op)
  }

  if(ALWAYS) {
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
  }

  if(ALWAYS) {
    #' case of zero cif
    cat("zero cif...", fill=TRUE)
    pmiH <- Everything(fitH, sparseOK=TRUE)
    pmiSH <- Everything(fitSH, sparseOK=TRUE)
    pmiHx <- Everything(fitHx, sparseOK=TRUE)
  }
  if(FULLTEST) {
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
  }
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
    cat(tag, fill=TRUE)
    #' compute all leverage+influence terms
    ppmInfluence(..., what="all", iScore=iS, iHessian=iH)
  }
  gogogo <- function(hdr, fit) {
    cat(hdr, fill=TRUE)
    force(fit)
    #' try all code options
    d <- gogo("a", fit)
    d <- gogo("b", fit, method="interpreted") 
    d <- gogo("c", fit, method="interpreted", entrywise=FALSE)
    d <- gogo("d", fit,                       entrywise=FALSE) 
    invisible(NULL)
  }
  gogogo("Offset model...",
         ippm(X ~ offset(ytoa), start=list(alpha=1), iterlim=40))
  gogogo("Offset model (logistic) ...", 
         ippm(X ~ offset(ytoa), start=list(alpha=1),
              method="logi", iterlim=40))
  gogogo("Offset+x model...", 
         ippm(X ~ x + offset(ytoa), start=list(alpha=1), iterlim=40))
  gogogo("Offset+x model (logistic) ...", 
         ippm(X ~ x + offset(ytoa), start=list(alpha=1),
              method="logi", iterlim=40))
  gogogo("Offset model Strauss ...", 
         ippm(X ~ offset(ytoa), Strauss(0.07), start=list(alpha=1), iterlim=40))
  gogogo("Offset model Strauss (logistic) ...", 
         ippm(X ~ offset(ytoa), Strauss(0.07), start=list(alpha=1),
              method="logi", iterlim=40))
  gogogo("Offset+x model Strauss ...", 
         ippm(X ~ x + offset(ytoa), Strauss(0.07), start=list(alpha=1),
              iterlim=40))
  gogogo("Offset+x model Strauss (logistic)...", 
         ippm(X ~ x + offset(ytoa), Strauss(0.07), start=list(alpha=1),
              method="logi", iterlim=40))
  #'
  set.seed(452)
  foo <- ppm(cells ~ 1, Strauss(0.15), method="ho", nsim=5)
  aa <- Everything(foo)

  #' Gradient and Hessian obtained by symbolic differentiation
  f <- deriv(expression((1+x)^a),
             "a", function.arg=c("x", "y", "a"),
             hessian=TRUE)
  #' check they can be extracted
  fit <- ippm(cells ~offset(f), start=list(a=0.7))
  Everything(fit)
})

reset.spatstat.options()
