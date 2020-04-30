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
#' spatstat/tests/package.R
#' Package information
#' $Revision: 1.2 $ $Date: 2020/04/30 05:23:52 $

if(FULLTEST) {
local({
  a <- bugfixes("book", show=FALSE)
  bugfixes(package="deldir")
})
}
## 
##    tests/percy.R
##
## Tests of Percus-Yevick approximations
##
##    $Revision: 1.3 $ $Date: 2020/04/30 05:23:52 $

if(FULLTEST) {
local({
  fit <- ppm(swedishpines ~1, DiggleGatesStibbard(6))
  K <- Kmodel(fit)
})
}

#'   tests/perspim.R
#'
#'   Check persp.im handling of NA, etc
#' 
#'   $Revision: 1.2 $  $Date: 2020/04/30 05:23:52 $

if(FULLTEST) {
local({
  set.seed(42)
  Z <- distmap(letterR, invert=TRUE)[letterR, drop=FALSE]
  X <- runifpoint(100, Frame(Z))
  M <- persp(Z, colin=Z, visible=TRUE, phi=50)
  perspPoints(X, Z=Z, M=M)
  P <- psp(c(2.360, 3.079, 2.211),
           c(0.934, 1.881, 2.184),
           c(2.337, 3.654, 3.274),
           c(1.829, 0.883, 2.093), window=letterR)
  perspSegments(P, Z=Z, M=M)
  
  persp(Z, colmap=rainbow)
  persp(Z, colmap=beachcolours, sealevel=mean(Z))
  persp(Z, colin=as.im(Z, dimyx=dim(Z)/4))
})
}
##
## tests/pixelgripes.R
##     Problems related to pixellation of windows
##
## $Revision: 1.5 $ $Date: 2020/04/30 05:23:52 $

if(FULLTEST) {
local({
  ## From Philipp Hunziker: bug in rNeymanScott (etc)
  ## Create an irregular window
  PM <- matrix(c(1,0,0.5,1,0,0), 3, 2, byrow=TRUE)
  P <- owin(poly=PM)
  ## Generate Matern points
  X <- rMatClust(50, 0.05, 5, win=P)
  ## Some distance function as a covariate
  distorigin <- function(x, y) { sqrt(x^2 + y^2) }
  ## No covariates: works fine
  fit0 <- kppm(X ~ 1, clusters="MatClust")
  Y0 <- simulate(fit0, retry=0)
  ## Covariates: Simulation fails
  fit1 <- kppm(X ~ distorigin, clusters="MatClust")
  Y1 <- simulate(fit1, retry=0)

  ## pixellate.ppp includes mapping from (x,y) to (row, col)
  Z <- pixellate(cells, savemap=TRUE)
  ind <- attr(Z, "map")
  m <- (as.matrix(Z))[ind]
  if(!all(m == 1)) stop("Coordinate mismatch in pixellate.ppp")
})
}
## 
## tests/polygons.R
##
##  $Revision: 1.5 $ $Date: 2020/04/30 05:23:52 $
##
if(ALWAYS) { # involves C code
local({
  co <- as.ppp(corners(letterR), letterR, check=FALSE)
  co[letterR]

  b <- letterR$bdry
  a <- sapply(b, xypolyselfint, yesorno=TRUE)
  a <- lapply(b, xypolyselfint, proper=TRUE)
  
  ## Simple example of self-crossing polygon
  x <- read.table("selfcross.txt", header=TRUE)
  y <- xypolyselfint(x)
})
}
