#
#        edgeTrans.R
#
#    $Revision: 1.15 $    $Date: 2016/04/25 02:34:40 $
#
#    Translation edge correction weights
#
#  edge.Trans(X)      compute translation correction weights
#                     for each pair of points from point pattern X 
#
#  edge.Trans(X, Y, W)   compute translation correction weights
#                        for all pairs of points X[i] and Y[j]
#                        (i.e. one point from X and one from Y)
#                        in window W
#
#  edge.Trans(X, Y, W, paired=TRUE)
#                        compute translation correction weights
#                        for each corresponding pair X[i], Y[i].
#
#  To estimate the K-function see the idiom in "Kest.R"
#
#######################################################################

edge.Trans <- function(X, Y=X, W=Window(X), exact=FALSE, paired=FALSE,
                       ..., 
                       trim=spatstat.options("maxedgewt"),
                       dx=NULL, dy=NULL,
                       give.rmax=FALSE,
                       gW = NULL) {
  given.dxdy <- !is.null(dx) && !is.null(dy)
  if(!given.dxdy) {
    ## dx, dy will be computed from X, Y
    X <- as.ppp(X, W)
    W <- X$window
    Y <- if(!missing(Y)) as.ppp(Y, W) else X
    nX <- X$n
    nY <- Y$n
    if(paired) {
      if(nX != nY)
        stop("X and Y should have equal length when paired=TRUE")
      dx <- Y$x - X$x
      dy <- Y$y - X$y
    } else {
      dx <- outer(X$x, Y$x, "-")
      dy <- outer(X$y, Y$y, "-")
    }
  } else {
    ## dx, dy given
    if(paired) {
      ## dx, dy are vectors
      check.nvector(dx)
      check.nvector(dy)
      stopifnot(length(dx) == length(dy))
    } else {
      ## dx, dy are matrices
      check.nmatrix(dx)
      check.nmatrix(dy)
      stopifnot(all(dim(dx) == dim(dy)))
      nX <- nrow(dx)
      nY <- ncol(dx)
    }
    stopifnot(is.owin(W))
  }
    
  ## For irregular polygons, exact evaluation is very slow;
  ## so use pixel approximation, unless exact=TRUE
  if(W$type == "polygonal" && !exact)
    W <- as.mask(W)

  ## compute
  if(!paired) {
    dx <- as.vector(dx)
    dy <- as.vector(dy)
  }
  switch(W$type,
         rectangle={
           ## Fast code for this case
           wide <- diff(W$xrange)
           high <- diff(W$yrange)
           weight <- wide * high / ((wide - abs(dx)) * (high - abs(dy)))
         },
         polygonal={
           ## This code is SLOW
           n <- length(dx)
           weight <- numeric(n)
           if(n > 0) {
             for(i in seq_len(n)) {
               Wshift <- shift(W, c(dx[i], dy[i]))
               weight[i] <- overlap.owin(W, Wshift)
             }
             weight <- area(W)/weight
           }
         },
         mask={
           ## compute set covariance of window
           if(is.null(gW)) gW <- setcov(W)
           ## evaluate set covariance at these vectors
           gvalues <- lookup.im(gW, dx, dy, naok=TRUE, strict=FALSE)
           weight <- area(W)/gvalues
         }
         )
  
  ## clip high values
  if(length(weight) > 0)
    weight <- pmin.int(weight, trim)

  if(!paired) 
    weight <- matrix(weight, nrow=nX, ncol=nY)

  if(give.rmax) 
    attr(weight, "rmax") <- rmax.Trans(W, gW)
  return(weight)
}

## maximum radius for translation correction
## = radius of largest circle centred at 0 contained in W + ^W

rmax.Trans <- function(W, g=setcov(W)) {
  ## calculate maximum permissible 'r' value
  ## for validity of translation correction
  W <- as.owin(W)
  if(is.rectangle(W)) 
    return(shortside(W))
  ## find support of set covariance
  if(is.null(g)) g <- setcov(W)
  eps <- 2 * max(1, max(g)) * .Machine$double.eps
  gsupport <- solutionset(g > eps)
  gboundary <- bdry.mask(gsupport)
  xy <- rasterxy.mask(gboundary, drop=TRUE)
  rmax <- with(xy, sqrt(min(x^2 + y^2)))
  return(rmax)
}

## maximum radius for rigid motion correction
## = radius of smallest circle centred at 0 containing W + ^W

rmax.Rigid <- function(X, g=setcov(Window(X))) {
  stopifnot(is.ppp(X) || is.owin(X))
  if(is.ppp(X))
    return(max(pairdist(X[chull(X)])))
  W <- X
  if(is.rectangle(W)) return(diameter(W))
  if(is.null(g)) g <- setcov(W)
  eps <- 2 * max(1, max(g)) * .Machine$double.eps
  gsupport <- solutionset(g > eps)
  gboundary <- bdry.mask(gsupport)
  xy <- rasterxy.mask(gboundary, drop=TRUE)
  rmax <- with(xy, sqrt(max(x^2 + y^2)))
  return(rmax)
}
