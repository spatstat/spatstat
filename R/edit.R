##   edit.R
##
##   Methods for 'edit'
##
##   $Revision: 1.3 $ $Date: 2015/04/19 06:14:21 $

edit.ppp <- local({

  edit.ppp <- function(name, ...) {
    X <- name
    df <- as.data.frame(X)
    df <- as.data.frame(lapply(df, as.num.or.char))
    Y <- edit(df, ...)
    Z <- as.ppp(Y, W=Window(X))
    return(Z)
  }

  as.num.or.char <- function(x) {
    if (is.character(x)) x else
    if (is.numeric(x)) {
      storage.mode(x) <- "double"
      x
    } else as.character(x)
  }

  edit.ppp
})

edit.im <- function(name, ...) {
  X <- name
  M <- transmat(as.matrix(X), from="spatstat", to="European")
  Y <- as.data.frame(M)
  Z <- edit(Y, ...)
  X[] <- transmat(as.matrix(Z), from="European", to="spatstat")
  return(X)
}

  
