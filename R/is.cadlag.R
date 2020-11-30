#'
#'   is.cadlag.R
#'
#'   Test whether a stepfun is cadlag/rcll
#'   (continue a droite; limites a gauche)
#'
#'   $Revision: 1.4 $ $Date: 2020/11/30 04:10:33 $


is.cadlag <- function (s) {
  stopifnot(is.stepfun(s))
  r <- knots(s)
  h <- s(r)
  n <- length(r)
  r1 <- c(r[-1L],r[n]+1)
  rm <- (r+r1)/2
  hm <- s(rm)
  isTRUE(all.equal(h,hm))
}

