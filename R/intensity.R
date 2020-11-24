#
# intensity.R
#
# Code related to intensity and intensity approximations
#
#  $Revision: 1.23 $ $Date: 2020/11/24 01:57:52 $
#

intensity <- function(X, ...) {
  UseMethod("intensity")
}

intensity.ppp <- function(X, ..., weights=NULL) {
  n <- npoints(X)
  a <- area(Window(X))
  if(is.null(weights)) {
    ## unweighted case - for efficiency
    if(is.multitype(X)) {
      mks <- marks(X)
      answer <- as.vector(table(mks))/a
      names(answer) <- levels(mks)
    } else answer <- n/a
    return(answer)
  }
  ## weighted case
  weights <- pointweights(X, weights=weights, parent=parent.frame())
  if(is.multitype(X)) {
    mks <- marks(X)
    answer <- as.vector(tapply(weights, mks, sum))/a
    answer[is.na(answer)] <- 0
    names(answer) <- levels(mks)
  } else {
    answer <- sum(weights)/a
  }
  return(answer)
}

intensity.splitppp <- function(X, ..., weights=NULL) {
  if(is.null(weights))
    return(sapply(X, intensity.ppp))
  if(is.expression(weights))
    return(sapply(X, intensity.ppp, weights=weights))
  if(is.numeric(weights)) {
    fsplit <- attr(X, "fsplit")
    n <- length(fsplit)
    check.nvector(weights, n)
    result <- mapply(intensity.ppp, X, weights=split(weights, fsplit))
    result <- simplify2array(result, higher=FALSE)
    return(result)
  }
  stop("Unrecognised format for weights")
}

