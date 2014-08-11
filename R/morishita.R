#
# morishita.R
#
#  $Revision: 1.5 $  $Date: 2013/04/25 06:37:43 $
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
  unitinfo <- summary(unitname(W))$axis
  do.call("plot.default",
          resolve.defaults(list(quadsize, mindex),
                           list(...),
                           list(xlim=c(0,a),
                                ylim=c(0,max(mindex)),
                                xlab=paste("Diameter of quadrat", unitinfo),
                                ylab="Morishita index",
                                main=paste("Morishita plot for", Xname))))
  abline(h=1, lty=2)
  return(invisible(NULL))
}
