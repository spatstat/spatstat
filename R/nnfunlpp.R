#
# nnfunlpp.R
#
#   method for 'nnfun' for class 'lpp'
#
#   $Revision: 1.2 $ $Date: 2012/10/18 10:55:16 $
#

nnfun.lpp <- local({

  nnfun.lpp <- function(X, ...) {
    stopifnot(inherits(X, "lpp"))
    force(X)
    L <- as.linnet(X)
    f <- function(x, y=NULL, seg=NULL, tp=NULL, ...) {
      # L is part of the environment
      Y <- as.lpp(x=x, y=y, seg=seg, tp=tp, L=L)
      i <- nncross.lpp(Y, X, what="which")
      return(i)
    }
    f <- linfun(f, L)
    attr(f, "explain") <- uitleggen
    return(f)
  }

  uitleggen <- function(x, ...) {
    cat("Nearest neighbour function for lpp object\n")
    X <-  get("X", envir=environment(x))
    print(X)
  }

  nnfun.lpp
})


