#'
#'  distcdf.R
#'
#' cdf of |X1-X2| when X1,X2 are iid uniform in W, etc
#'
#'  $Revision: 1.10 $  $Date: 2016/02/11 10:17:12 $
#'

distcdf <- function(W, V=W, ..., dW=1, dV=dW, nr=1024, regularise=TRUE) {
  reflexive <- missing(V) && missing(dV)
  diffuse <- is.owin(W) && is.owin(V)
  uniformW <- identical(dW, 1)
  uniformV <- identical(dV, 1)
  uniform <- uniformW && uniformV

  if(is.owin(W)) {
    W <- as.mask(as.owin(W), ...)
    dW <- as.im(dW, W=W)
  } else if(is.ppp(W)) {
    if(uniformW) {
      #' discrete uniform distribution on W
      dW <- pixellate(W, ...)
    } else {
      #' dW should be a weight or vector of weights
      if(!is.vector(dW) || !is.numeric(dW))
        stop("If W is a point pattern, dW should be a vector of weights")
      if(length(dW) == 1L) {
        dW <- rep(dW, npoints(W))
      } else stopifnot(length(dW) == npoints(W))
      dW <- pixellate(W, weights=dW, ...)
    }
  } else stop("W should be a point pattern or a window")
  
  if(is.owin(V)) {
    V <- as.mask(as.owin(V), ...)
    dV <- as.im(dV, W=V)
  } else if(is.ppp(V)) {
    if(uniformV) {
      #' discrete uniform distribution on V
      dV <- pixellate(V, ...)
    } else {
      #' dV should be a weight or vector of weights
      if(!is.vector(dV) || !is.numeric(dV))
        stop("If V is a point pattern, dV should be a vector of weights")
      if(length(dV) == 1L) {
        dV <- rep(dV, npoints(V))
      } else stopifnot(length(dV) == npoints(V))
      dV <- pixellate(V, weights=dV, ...)
    }
  } else stop("V should be a point pattern or a window")

  if(!uniformW && min(dW) < 0) 
    stop("Negative values encountered in dW")
  
  if(!uniformV && min(dV) < 0) 
    stop("Negative values encountered in dV")

  #' compute
  if(diffuse && uniform) {
    #' uniform distributions on windows 
    g <- if(reflexive) setcov(W, ...) else setcov(W, V, ...)
  } else {
    g <- if(reflexive) imcov(dW) else imcov(dW, dV)
  }
  r <- as.im(function(x,y) { sqrt(x^2 + y^2) }, g)
  rvals <- as.vector(as.matrix(r))
  gvals <- as.vector(as.matrix(g))
  rgrid <- seq(0, max(rvals), length=nr)
  dr <- max(rvals)/(nr-1)
  h <- whist(rvals, breaks=rgrid, weights=gvals/sum(gvals))
  ch <- c(0,cumsum(h))
  #' regularise at very short distances
  if(regularise) {
    sevenpix <- 7 * with(r, max(xstep, ystep))
    ii <- round(sevenpix/dr)
    ch[1:ii] <- ch[ii] * (rgrid[1:ii]/rgrid[ii])^2
  }
  #' ok
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
  do.call(plot,
          resolve.defaults(list(g),
                             list(...),
                             list(main=xname)))
  abline(v=ropt, lty=3)
  abline(h=f, lty=3)
  invisible(NULL)
}

