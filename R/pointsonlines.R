#
#   pointsonlines.R
#
# place points at regular intervals along line segments
#
#   $Revision: 1.9 $  $Date: 2020/03/16 10:28:51 $
#

pointsOnLines <- function(X, eps=NULL, np=1000, shortok=TRUE) {
  stopifnot(is.psp(X))
  len <- lengths_psp(X)
  nseg <- length(len)
  if(is.null(eps)) {
    stopifnot(is.numeric(np) && length(np) == 1)
    stopifnot(is.finite(np) && np > 0)
    eps <- sum(len)/np
  } else {
    stopifnot(is.numeric(eps) && length(eps) == 1)
    stopifnot(is.finite(eps) && eps > 0)
  }
  # initialise
  Xdf    <- as.data.frame(X)
  xmid <- with(Xdf, (x0+x1)/2)
  ymid <- with(Xdf, (y0+y1)/2)
  # handle very short segments
#  allsegs <- 1:nseg
  if(any(short <- (len <= eps)) && shortok) {
    # very short segments: use midpoints
    Z <- data.frame(x = xmid[short], y = ymid[short], seg=which(short), tp=0.5)
  } else Z <- data.frame(x=numeric(0), y=numeric(0),
                         seg=integer(0), tp=numeric(0))
  # handle other segments
  for(i in (1:nseg)[!short]) {
    # divide segment into pieces of length eps
    # with shorter bits at each end
    leni <- len[i]
    nwhole <- floor(leni/eps)
    if(leni/eps - nwhole < 0.5 && nwhole > 2)
      nwhole <- nwhole - 1
    rump <- (leni - nwhole * eps)/2
    brks <- c(0, rump + (0:nwhole) * eps, leni)
    nbrks <- length(brks)
    # points at middle of each piece
    ss <- (brks[-1] + brks[-nbrks])/2
    tp <- ss/leni
    x <- with(Xdf, x0[i] + tp * (x1[i]-x0[i]))
    y <- with(Xdf, y0[i] + tp * (y1[i]-y0[i]))
    Z <- rbind(Z, data.frame(x=x, y=y, seg=i, tp=tp))
  }
  result <- as.ppp(Z[,c("x","y")], W=X$window)
  attr(result, "map") <- Z[,c("seg", "tp")]
  return(result)
}
