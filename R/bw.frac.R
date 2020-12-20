#'
#'    bw.frac.R
#'
#'   $Revision: 1.2 $  $Date: 2020/12/19 05:25:06 $

bw.frac <- function(X, ..., f=1/4) {
  X <- as.owin(X)
  g <- distcdf(X, ...)
  r <- with(g, .x)
  Fr <- with(g, .y)
  iopt <- min(which(Fr >= f))
  ropt <- r[iopt]
  attr(ropt, "f") <- f
  attr(ropt, "g") <- g
  class(ropt) <- c("bw.frac", class(ropt))
  return(ropt)
}

print.bw.frac <- function(x, ...) {
  print(as.numeric(x), ...)
}

plot.bw.frac <- function(x, ...) {
  xname <- short.deparse(substitute(x))
  g <- attr(x, "g")
  f <- attr(x, "f")
  ropt <- as.numeric(x)
  dont.complain.about(g)
  do.call(plot,
          resolve.defaults(list(quote(g)),
                             list(...),
                             list(main=xname)))
  abline(v=ropt, lty=3)
  abline(h=f, lty=3)
  invisible(NULL)
}

