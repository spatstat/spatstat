fourierbasis <- function(x, k, win = boxx(rep(list(0:1), ncol(k)))) {
  x <- as.matrix(x)
  k <- as.matrix(k)
  if (nrow(k) == 0 | nrow(x) == 0) 
    return(complex())
  d <- ncol(x)
  if (ncol(k) != d) 
    stop("Arguments x and k must have the same number of columns.")
  win <- as.boxx(win)
  boxlengths <- as.numeric(win$ranges[2L, ] - win$ranges[1L, ])
  if (length(boxlengths) != d) 
    stop("The box dimension differs from the number of columns in x and k")
  rslt <- exp(2 * pi * (0+1i) * outer(k[, 1L], x[, 1L]/boxlengths[1L]))
  if (d > 1) {
    for (i in 2:d) {
      rslt <- rslt * exp(2 * pi * (0+1i) * outer(k[, i], x[, i]/boxlengths[i]))
    }
  }
  return(rslt/prod(boxlengths))
}
