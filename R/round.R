#
#   round.R
#
#   discretisation of coordinates
#
#   $Revision: 1.5 $  $Date: 2013/01/09 03:13:10 $

round.ppp <- round.pp3 <- round.ppx <- function(x, digits=0) {
  coords(x) <- round(as.matrix(coords(x)), digits=digits)
  return(x)
}

rounding <- function(x) {
  UseMethod("rounding")
}

rounding.ppp <- rounding.pp3 <- rounding.ppx <- function(x) {
  rounding(as.matrix(coords(x)))
}

rounding.default <- function(x) {
  # works for numeric, complex, matrix etc
  if(all(x == 0))
    return(NULL)
  if(identical(all.equal(x, round(x)), TRUE)) { 
    # integers: go up
    k <- 0
    smallk <- -log10(.Machine$double.xmax)
    repeat {
      if(k < smallk || !identical(all.equal(x, round(x, k-1)), TRUE))
        return(k)
      k <- k-1
    }
  } else {
    # not integers: go down
    k <- 1
    bigk <- -log10(.Machine$double.eps)
    repeat {
      if(k > bigk || identical(all.equal(x, round(x, k)), TRUE))
        return(k)
      k <- k+1
    }
  }
}
