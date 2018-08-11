#'
#'     headtail.R
#'
#'   Methods for head() and tail()
#'
#'   $Revision: 1.1 $  $Date: 2016/12/20 01:11:29 $

head.tess <- head.psp <- head.ppx <- head.ppp <- function(x, n=6L, ...) {
  stopifnot(length(n) == 1L)
  xlen <- nobjects(x)
  n <- if (n < 0L) max(xlen + n, 0L) else min(n, xlen)
  x[seq_len(n)]
}

tail.tess <- tail.psp <- tail.ppx <- tail.ppp <- function (x, n = 6L, ...) {
  stopifnot(length(n) == 1L)
  xlen <- nobjects(x)
  n <- if (n < 0L) max(xlen + n, 0L) else min(n, xlen)
  x[seq.int(to = xlen, length.out = n)]
}


