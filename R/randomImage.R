#'
#' randomImage.R
#'
#' Functions for generating random images
#' 
#'    $Revision: 1.1 $  $Date: 2015/03/23 10:44:04 $
#'
#'

rnoise <- function(rgen=runif, w=square(1), ...) {
  a <- do.call.matched(as.mask, list(w=w, ...), sieve=TRUE)
  W <- a$result
  argh <- a$otherargs
  Z <- as.im(W)
  n <- sum(W$m)
  Z[] <- do.call(rgen, append(list(n=n), argh))
  return(Z)
}

  
  
