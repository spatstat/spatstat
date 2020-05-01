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
#'    tests/quadschemes.R
#'
#'    Quadrature schemes, dummy points etc
#' 
#'   $Revision: 1.7 $ $Date: 2020/05/01 02:42:58 $
#'

if(FULLTEST) {
local({
  ##  class 'quad' 
  qu <- quadscheme(cells)
  qm <- quadscheme(amacrine)
  plot(qu)
  plot(qm)
  is.multitype(qu)
  is.multitype(qm)
  a <- param.quad(qu)
  a <- param.quad(qm)
  a <- equals.quad(qu)
  a <- equals.quad(qm)
  a <- domain(qu)
  unitname(qu) <- c("Furlong", "Furlongs")
  
  ## utilities
  b <- cellmiddles(square(1), 3, 4)
  b <- cellmiddles(letterR, 3, 4, distances=FALSE)
  b <- cellmiddles(letterR, 3, 4, distances=TRUE)
  v <- tilecentroids(square(1), 3, 4)
  v <- tilecentroids(letterR, 3, 4)
  n <- default.n.tiling(cells)
  n <- default.n.tiling(cells, nd=4)
  n <- default.n.tiling(cells, ntile=4)
  n <- default.n.tiling(cells, ntile=4, quasi=TRUE)

  ## quadrature weights - special cases
  X <- runifpoint(10, as.mask(letterR))
  gr <- gridweights(X, ntile=12, npix=7) # causes warnings about zero digital area
  
  ## plot.quad 
  plot(quadscheme(cells, method="dirichlet", nd=7),              tiles=TRUE)
  plot(quadscheme(cells, method="dirichlet", nd=7, exact=FALSE), tiles=TRUE)
  
  ## logistic
  d <- quadscheme.logi(cells, logi.dummy(cells, "binomial"))
  print(summary(d))
  d <- quadscheme.logi(cells, logi.dummy(cells, "poisson"))
  print(summary(d))
  d <- quadscheme.logi(cells, logi.dummy(cells, "grid"))
  print(summary(d))
  d <- quadscheme.logi(cells, logi.dummy(cells, "transgrid"))
  print(summary(d))
  d <- quadscheme.logi(amacrine,
                       logi.dummy(amacrine, "binomial", mark.repeat=TRUE))
  print(summary(d))
  d <- quadscheme.logi(amacrine,
                       logi.dummy(amacrine, "poisson", mark.repeat=FALSE))
  print(summary(d))
})
}
