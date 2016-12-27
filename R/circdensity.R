#'
#'   circdensity.R
#'
#' Kernel smoothing for circular data
#'
#'   $Revision: 1.3 $ $Date: 2014/12/04 06:49:20 $

circdensity <- function(x, sigma="nrd0", ..., bw=NULL,
                        weights=NULL,
                        unit=c("degree", "radian")) {
  xname <- short.deparse(substitute(x))
  missu <- missing(unit)
  if(missing(sigma) && !is.null(bw))
    sigma <- bw
  unit <- match.arg(unit)
  unit <- validate.angles(x, unit, missu)
  FullCircle <- switch(unit, degree = 360, radian = 2*pi)
  if(is.character(sigma)) {
    sigma <- switch(sigma,
                     bcv  = bw.bcv,
                     nrd  = bw.nrd,
                     nrd0 = bw.nrd0,
                     SJ   = bw.SJ,
                     ucv  = bw.ucv,
                     get(paste0("bw.", sigma), mode="function"))
  }
  if(is.function(sigma)) {
    sigma <- sigma(x)
    if(!(is.numeric(sigma) && length(sigma) == 1L && sigma > 0))
      stop("Bandwidth selector should return a single positive number")
  }
  check.1.real(sigma)
  #' replicate data
  x <- x %% FullCircle
  xx <- c(x - FullCircle, x, x + FullCircle)
  #' replicate weights
  if(!is.null(weights)) {
    stopifnot(length(weights) == length(x))
    weights <- rep(weights, 3)/3
  }
  #' smooth
  z <- do.call(density.default,
               resolve.defaults(list(x=xx, bw=sigma, weights=weights),
                                list(...),
                                list(from=0, to=FullCircle)))
  z$y <- 3 * z$y
  z$data.name <- xname
  return(z)
}
