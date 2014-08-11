#
#        edgeTrans.R
#
#    $Revision: 1.11 $    $Date: 2011/05/18 01:51:52 $
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
#  To estimate the K-function see the idiom in "Kest.S"
#
#######################################################################

edge.Trans <- function(X, Y=X, W=X$window, exact=FALSE, paired=FALSE,
                       trim=spatstat.options("maxedgewt")) {

  X <- as.ppp(X, W)

  W <- X$window
  x <- X$x
  y <- X$y
  nX <- X$n
  
  Y <- as.ppp(Y, W)
  xx <- Y$x
  yy <- Y$y
  nY <- Y$n
  
  if(paired && (nX != nY))
    stop("X and Y should have equal length when paired=TRUE")
  
  # For irregular polygons, exact evaluation is very slow;
  # so use pixel approximation, unless exact=TRUE
  if(W$type == "polygonal" && !exact)
    W <- as.mask(W)

  switch(W$type,
         rectangle={
           # Fast code for this case
           wide <- diff(W$xrange)
           high <- diff(W$yrange)
           if(!paired) {
             DX <- abs(outer(x,xx,"-"))
             DY <- abs(outer(y,yy,"-"))
           } else {
             DX <- abs(xx - x)
             DY <- abs(yy - y)
           }
           weight <- wide * high / ((wide - DX) * (high - DY))
         },
         polygonal={
           # This code is SLOW
           a <- area.owin(W)
           if(!paired) {
             weight <- matrix(, nrow=nX, ncol=nY)
             if(nX > 0 && nY > 0) {
               for(i in seq_len(nX)) {
                 X.i <- c(x[i], y[i])
                 for(j in seq_len(nY)) {
                   shiftvector <- X.i - c(xx[j],yy[j])
                   Wshift <- shift(W, shiftvector)
                   b <- overlap.owin(W, Wshift)
                   weight[i,j] <- a/b
                 }
               }
             }
           } else {
             weight <- numeric(nX)
             if(nX > 0) {
               for(i in seq_len(nX)) {
                 shiftvector <- c(x[i],y[i]) - c(xx[i],yy[i])
                 Wshift <- shift(W, shiftvector)
                 b <- overlap.owin(W, Wshift)
                 weight[i] <- a/b
               }
             }
           }
         },
         mask={
           # make difference vectors
           if(!paired) {
             DX <- outer(x,xx,"-")
             DY <- outer(y,yy,"-")
           } else {
             DX <- x - xx
             DY <- y - yy
           }
           # compute set covariance of window
           g <- setcov(W)
           # evaluate set covariance at these vectors
           gvalues <- lookup.im(g, as.vector(DX), as.vector(DY),
                                naok=TRUE, strict=FALSE)
           if(!paired) 
             # reshape
             gvalues <- matrix(gvalues, nrow=nX, ncol=nY)
           weight <- area.owin(W)/gvalues
         }
         )
  # clip high values
  if(length(weight) > 0)
    weight <- pmin.int(weight, trim)
  if(!paired) 
    weight <- matrix(weight, nrow=nX, ncol=nY)
  return(weight)
}
