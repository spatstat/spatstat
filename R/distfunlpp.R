#
# distfunlpp.R
#
#   method for 'distfun' for class 'lpp'
#
#   $Revision: 1.3 $ $Date: 2018/04/23 04:52:17 $
#

distfun.lpp <- local({
  
  distfun.lpp <- function(X, ..., k=1) {
    stopifnot(inherits(X, "lpp"))
    force(X)
    force(k)
    stopifnot(length(k) == 1)
    L <- as.linnet(X)
    f <- function(x, y=NULL, seg=NULL, tp=NULL, ...) {
      # L is part of the environment
      Y <- as.lpp(x=x, y=y, seg=seg, tp=tp, L=L)
      d <- nncross.lpp(Y, X, what="dist", k=k)
      return(d)
    }
    f <- linfun(f, L)
    assign("k", k, envir=environment(f))
    assign("X", X, envir=environment(f))
    attr(f, "explain") <- uitleggen
    attr(f, "extrargs") <- list(k=k)
    return(f)
  }

  uitleggen <- function(x, ...) {
    splat("Distance function for lpp object")
    envx <- environment(x)
    k <-  get("k", envir=envx)
    if(k != 1L)
      splat("Yields distance to", ordinal(k), "nearest point")
    X <-  get("X", envir=envx)
    print(X)
  }

  distfun.lpp
})




