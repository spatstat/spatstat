#
# nnfunlpp.R
#
#   method for 'nnfun' for class 'lpp'
#
#   $Revision: 1.3 $ $Date: 2019/09/16 10:14:18 $
#

nnfun.lpp <- local({

  nnfun.lpp <- function(X, ..., k=1, value=c("index", "mark")) {
    stopifnot(inherits(X, "lpp"))
    force(X)
    force(k)
    value <- match.arg(value)
    L <- as.linnet(X)
    switch(value,
           index = {
             fi <- function(x, y=NULL, seg=NULL, tp=NULL, ...) {
               ## L is part of the environment
               Y <- as.lpp(x=x, y=y, seg=seg, tp=tp, L=L)
               i <- nncross.lpp(Y, X, what="which", k=k)
               return(i)
             }
             f <- linfun(fi, L)
           },
           mark = {
             stopifnot(is.marked(X))
             marx <- as.data.frame(marks(X))[,1]
             fm <- function(x, y=NULL, seg=NULL, tp=NULL, ...) {
               Y <- as.lpp(x=x, y=y, seg=seg, tp=tp, L=L)
               i <- nncross.lpp(Y, X, what="which", k=k)
               return(marx[i])
             }
             f <- linfun(fm, L)
           })
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
    v <- mget("value", envir=env, ifnotfound=list(NULL))[[1L]]
    splat("Function returns the",
          if(identical(v, "mark")) "mark value" else "index",
          "of the neighbour")
  }

  nnfun.lpp
})


