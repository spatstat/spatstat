#
#  distcdf.R
#
# cdf of |X1-X2| when X1,X2 are iid uniform in W, etc
#
#  $Revision: 1.2 $  $Date: 2012/10/29 10:00:07 $
#

distcdf <- function(W, V=W, ..., dW=1, dV=dW, nr=1024) {
  if(missing(dW) && missing(dV)) {
    # uniform distributions
    g <- if(missing(V)) setcov(W, ...) else setcov(W, V, ...)
  } else {
    # nonuniform distributions
    W <- as.mask(as.owin(W), ...)
    V <- if(!missing(V)) as.mask(as.owin(V), ...) else W
    dW <- as.im(dW, W=W)
    if(!missing(V) || !missing(dV))
      dV <- as.im(dV, W=V)
    g <- imcov(dW, dV)
  }
  r <- as.im(function(x,y) { sqrt(x^2 + y^2) }, g)
  rvals <- as.vector(as.matrix(r))
  gvals <- as.vector(as.matrix(g))
  rgrid <- seq(0, max(rvals), length=nr)
  h <- whist(rvals, breaks=rgrid, weights=gvals/sum(gvals))
  ch <- c(0,cumsum(h))
  result <- fv(data.frame(r=rgrid, f=ch),
                "r", quote(CDF(r)),
               "f", , range(rvals), c("r","%s(r)"),
               c("Interpoint distance","Cumulative probability"),
               fname="CDF")
  return(result)
}

bw.frac <- function(X, ..., f=1/4) {
  g <- distcdf(X, ...)
  r <- with(g, .x)
  Fr <- with(g, .y)
  iopt <- min(which(Fr >= f))
  ropt <- r[iopt]
  attr(ropt, "f") <- f
  attr(ropt, "g") <- g
  class(ropt) <- c("bw.frac", class(ropt))
  return(ropt)
}

print.bw.frac <- function(x, ...) {
  print(as.numeric(x), ...)
}

plot.bw.frac <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  g <- attr(x, "g")
  f <- attr(x, "f")
  ropt <- as.numeric(x)
  do.call("plot",
          resolve.defaults(list(g),
                             list(...),
                             list(main=xname)))
  abline(v=ropt, lty=3)
  abline(h=f, lty=3)
  invisible(NULL)
}

