#
# morishita.R
#
#  $Revision: 1.6 $  $Date: 2014/06/03 02:58:09 $
#

miplot <- function(X, ...) {
  Xname <- short.deparse(substitute(X))
  X <- as.ppp(X)
  W <- X$window
  N <- X$n
  if(W$type != "rectangle")
    stop("Window of X is not a rectangle - Morishita index undefined")
  a <- min(diff(W$xrange), diff(W$yrange))
  maxnquad <- floor(a/mean(nndist(X)))
  if(maxnquad <= 1)
    stop("Not enough points for a Morishita plot")
  mindex <- numeric(maxnquad)
  for(nquad in 1:maxnquad) {
    qq <- quadratcount(X, nquad, nquad)
    tt <- as.vector(as.table(qq))
    mindex[nquad] <- length(tt) * sum(tt * (tt-1))/(N*(N-1))
  }

  quadsize <- diameter(W)/(1:maxnquad)
  ok <- (quadsize <= a)
  quadsize <- quadsize[ok]
  mindex   <- mindex[ok]
  
  unitinfo <- summary(unitname(W))$axis
  do.call("plot.default",
          resolve.defaults(list(quadsize, mindex),
                           list(...),
                           list(xlim=c(0,max(quadsize)),
                                ylim=c(0,max(1, mindex)),
                                xlab=paste("Diameter of quadrat", unitinfo),
                                ylab="Morishita index",
                                main=paste("Morishita plot for", Xname))))
  abline(h=1, lty=2)
  return(invisible(NULL))
}
