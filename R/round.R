#
#   round.R
#
#   discretisation of coordinates
#
#   $Revision: 1.6 $  $Date: 2019/02/20 03:34:50 $

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
  if(isTRUE(all.equal(x, round(x)))) { 
    # integers: go up
    k <- 0
    smallk <- -log10(.Machine$double.xmax)
    repeat {
      if(k < smallk || !isTRUE(all.equal(x, round(x, k-1))))
        return(k)
      k <- k-1
    }
  } else {
    # not integers: go down
    k <- 1
    bigk <- -log10(.Machine$double.eps)
    repeat {
      if(k > bigk || isTRUE(all.equal(x, round(x, k))))
        return(k)
      k <- k+1
    }
  }
}
