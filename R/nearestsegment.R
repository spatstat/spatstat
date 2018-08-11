#
#  nearestsegment.R
#
#  $Revision: 1.12 $  $Date: 2018/03/07 01:56:36 $
#
# Given a point pattern X and a line segment pattern Y,
# for each point x of X, determine which segment of Y is closest to x
# and find the point on Y closest to x.
#

nearestsegment <- function(X,Y) {
  return(ppllengine(X,Y,"identify"))
}

project2segment <- function(X, Y) {
  return(ppllengine(X,Y,"project"))
}
  
ppllengine <- function(X, Y, action="project", check=FALSE) {
  stopifnot(is.ppp(X))
  stopifnot(is.psp(Y))
  stopifnot(action %in% c("distance", "identify", "project"))
  # deal with empty patterns
  if(X$n == 0) {
    nowt <- numeric(0)
    none <- integer(0)
    switch(action,
           identify = return(none),
           distance = return(list(dist=nowt, which=none)),
           project  = return(list(Xproj=X, mapXY=none, d=nowt, tp=nowt)))
  }
  if(Y$n == 0)
    stop("Segment pattern Y contains 0 segments; projection undefined")
  #              
  XX <- as.matrix(as.data.frame(unmark(X)))
  YY <- as.matrix(as.data.frame(unmark(Y)))
  # determine which segment lies closest to each point
  huge <- max(diameter(as.rectangle(as.owin(X))),
              diameter(as.rectangle(as.owin(Y))))
  d <- distppllmin(XX, YY, huge^2)
  mapXY <- d$min.which
  if(action == "identify")
    return(mapXY)
  else if(action == "distance") 
    return(data.frame(dist=d$min.d, which=mapXY))
  
  # combine relevant rows of data
  alldata <- as.data.frame(cbind(XX, YY[mapXY, ,drop=FALSE]))
  colnames(alldata) <- c("x", "y", "x0", "y0", "x1", "y1")
  # coordinate geometry
  dx <- with(alldata, x1-x0)
  dy <- with(alldata, y1-y0)
  leng <- sqrt(dx^2 + dy^2)
  # rotation sines & cosines (may include 0/0)
  co <- dx/leng
  si <- dy/leng
  # vector to point from first endpoint of segment
  xv <- with(alldata, x - x0)
  yv <- with(alldata, y - y0)
  # rotate coordinate system so that x axis is parallel to line segment
  xpr <- xv * co + yv * si
#  ypr <- - xv * si + yv * co
  # determine whether projection is an endpoint or interior point of segment
  ok <- is.finite(xpr)
  left <- !ok | (xpr <= 0)
  right <- ok &  (xpr >= leng)
  # location of projected point in rotated coordinates
  xr <- with(alldata, ifelseAX(left, 0, ifelseXY(right, leng, xpr)))
  # back to standard coordinates
  xproj <- with(alldata, x0 + ifelseXB(ok, xr * co, 0))
  yproj <- with(alldata, y0 + ifelseXB(ok, xr * si, 0))
  Xproj <- ppp(xproj, yproj, window=X$window, marks=X$marks, check=check)
  # parametric coordinates
  tp <- xr/leng
  tp[!is.finite(tp)] <- 0
  # 
  return(list(Xproj=Xproj, mapXY=mapXY, d=d$min.d, tp=tp))
}

