#
# nnfunlpp.R
#
#   method for 'nnfun' for class 'lpp'
#
#   $Revision: 1.2 $ $Date: 2016/08/21 04:33:47 $
#

nnfun.lpp <- local({

  nnfun.lpp <- function(X, ..., k=1) {
    stopifnot(inherits(X, "lpp"))
    force(X)
    force(k)
    L <- as.linnet(X)
    f <- function(x, y=NULL, seg=NULL, tp=NULL, ...) {
      # L is part of the environment
      Y <- as.lpp(x=x, y=y, seg=seg, tp=tp, L=L)
      i <- nncross.lpp(Y, X, what="which", k=k)
      return(i)
    }
    f <- linfun(f, L)
    attr(f, "explain") <- uitleggen
    return(f)
  }

  uitleggen <- function(x, ...) {
    env <- environment(attr(x, "f"))
    X <- get("X", envir=env)
    k <- get("k", envir=env)
    if(identical(k, 1)) {
      cat("Nearest-neighbour function for lpp object\n")
    } else {
      cat("k-th nearest neighbour function for lpp object\n")
      cat(paste("k =", commasep(k), "\n"))
    }
    print(X)
  }

  nnfun.lpp
})


