#
# distfunlpp.R
#
#   method for 'distfun' for class 'lpp'
#
#   $Revision: 1.4 $ $Date: 2012/10/21 02:52:11 $
#

distfun.lpp <- local({
  
  distfun.lpp <- function(X, ...) {
    stopifnot(inherits(X, "lpp"))
    force(X)
    L <- as.linnet(X)
    f <- function(x, y=NULL, seg=NULL, tp=NULL, ...) {
      # L is part of the environment
      Y <- as.lpp(x=x, y=y, seg=seg, tp=tp, L=L)
      d <- nncross.lpp(Y, X, what="dist")
      return(d)
    }
    f <- linfun(f, L)
    attr(f, "explain") <- uitleggen
    return(f)
  }

  uitleggen <- function(x, ...) {
    cat("Distance function for lpp object\n")
    X <-  get("X", envir=environment(x))
    print(X)
  }

  distfun.lpp
})




