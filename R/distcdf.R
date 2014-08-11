#
#  distcdf.R
#
# cdf of |X1-X2| when X1,X2 are iid uniform in W, etc
#
#  $Revision: 1.6 $  $Date: 2013/08/26 09:51:51 $
#

distcdf <- function(W, V=W, ..., dW=1, dV=dW, nr=1024) {
  V.given <- !missing(V)
  dV.given <- !missing(dV)
  reflexive <- !V.given && !dV.given
  uniform <- missing(dW) && missing(dV)
  Wdud <- is.ppp(W) && missing(dW)
  Vdud <- V.given && is.ppp(V) && missing(dV) && missing(dW)
  diffuse <- !Wdud && !Vdud
  force(dV)
  # handle case where W or V is a point pattern
  # interpreted as a discrete uniform distribution
  if(Wdud) dW <- pixellate(W, ...)
  if(Vdud) dV <- pixellate(V, ...)
  # 
  if(diffuse && uniform) {
    # uniform distributions on windows 
    g <- if(reflexive) setcov(W, ...) else setcov(W, V, ...)
  } else {
    # nonuniform distribution(s)
    if(!Wdud) {
      W <- as.mask(as.owin(W), ...)
      dW <- as.im(dW, W=W)
    }
    if(reflexive) {
      g <- imcov(dW)
    } else {
      if(!Vdud) {
        V <- if(V.given) as.mask(as.owin(V), ...) else W
        dV <- as.im(dV, W=V)
      }
      g <- imcov(dW, dV)
    }
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
  X <- as.owin(X)
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

