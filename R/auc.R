##
## auc.R
##
##  Calculate area under ROC curve
##
## $Revision: 1.1 $ $Date: 2015/04/13 09:31:21 $

auc <- function(model) {
  if(is.multitype(model))
    stop("Sorry, not yet implemented for multitype models")
  if(is.stationary(model)) {
    aobs <- atheo <- 1/2
  } else {
    lambda <- intensity(model)
    Fl <- ecdf(lambda[])
    lambda <- as.im(lambda, Window(model))
    X <- data.ppm(model)
    lamX <- lambda[X]
    aobs <- mean(Fl(lamX))
    atheo <- mean(lambda[] * Fl(lambda[]))/mean(lambda)
  }
  result <- c(aobs, atheo)
  names(result) <- c("obs", "theo")
  return(result)
}

