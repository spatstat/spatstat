#' spatstat/tests/package.R
#' Package information
#' $Revision$ $Date$

require(spatstat)
local({
  a <- bugfixes("book", show=FALSE)
  bugfixes(package="deldir")
})
## 
##    tests/percy.R
##
## Tests of Percus-Yevick approximations
##
##    $Revision: 1.2 $ $Date: 2015/12/29 08:54:49 $

require(spatstat)
local({
  fit <- ppm(swedishpines ~1, DiggleGatesStibbard(6))
  K <- Kmodel(fit)
})

#'   tests/perspim.R
#'
#'   Check persp.im handling of NA, etc
#' 
#'   $Revision: 1.1 $  $Date: 2016/08/27 02:53:35 $

require(spatstat)

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
##
## tests/pixelgripes.R
##     Problems related to pixellation of windows
##
## $Revision: 1.4 $ $Date: 2018/10/10 08:04:10 $

require(spatstat)
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
})

local({
  ## pixellate.ppp includes mapping from (x,y) to (row, col)
  Z <- pixellate(cells, savemap=TRUE)
  ind <- attr(Z, "map")
  m <- (as.matrix(Z))[ind]
  if(!all(m == 1)) stop("Coordinate mismatch in pixellate.ppp")
})

## 
## tests/polygons.R
##
##  $Revision: 1.4 $ $Date: 2019/11/03 03:18:05 $
##
require(spatstat)
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

