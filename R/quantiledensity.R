#'
#'   quantiledensity.R
#'
#'  quantile method for class 'density'
#'
#'  Also a CDF from a 'density'
#' 
#'  $Revision: 1.3 $ $Date: 2015/09/01 11:53:15 $

quantile.density <- local({

  quantile.density <- function(x, probs = seq(0, 1, 0.25), names = TRUE, ...,
                               warn=TRUE) {
    stopifnot(inherits(x, "density"))
    #' check whether density estimate was restricted to an interval
    if(warn && is.call(cl <- x$call) && any(c("from", "to") %in% names(cl)))
      warning(paste("Density was normalised within the computed range",
                    "of x values", prange(c(cl$from, cl$to))),
              call.=FALSE)
    #' validate probs
    eps <- 100 * .Machine$double.eps
    if(any((p.ok <- !is.na(probs)) & (probs < -eps | probs > 1 + eps))) 
      stop("'probs' outside [0,1]")
    if (na.p <- any(!p.ok)) {
      o.pr <- probs
      probs <- probs[p.ok]
      probs <- pmax(0, pmin(1, probs))
    }
    np <- length(probs)
    qs <- rep(NA_real_, np)
    if (np > 0) {
      #' extract density values 
      xx <- x$x
      yy <- x$y
      nn <- length(xx)
      #' integrate, normalise
      Fx <- cumsum(yy * c(0, diff(xx)))
      Fx <- Fx/Fx[nn]
      #' quantile
      for(j in 1:np) {
        ii <- min(which(Fx >= probs[j]))
        if(!is.na(ii) && ii >= 1 && ii <= nn) 
          qs[j] <- xx[ii]
      }
      if (names && np > 0L) {
        names(qs) <- format_perc(probs)
      }
    }
    if (na.p) {
      o.pr[p.ok] <- qs
      names(o.pr) <- rep("", length(o.pr))
      names(o.pr)[p.ok] <- names(qs)
      return(o.pr)
    } else return(qs)
  }

  format_perc <- function (x, digits = max(2L, getOption("digits")),
                           probability = TRUE, use.fC = length(x) < 100, ...) {
    if (length(x)) {
      if (probability) x <- 100 * x
      paste0(if (use.fC) 
             formatC(x, format = "fg", width = 1, digits = digits)
      else format(x, trim = TRUE, digits = digits, ...), "%")
    }
    else character(0)
  }

  quantile.density
})


CDF <- function(f, ...) {
  UseMethod("CDF")
}

CDF.density <- function(f, ..., warn=TRUE) {
  stopifnot(inherits(f, "density"))
  #' check whether density estimate was restricted to an interval
  if(warn && is.call(cl <- f$call) && any(c("from", "to") %in% names(cl)))
    warning(paste("Density was normalised within the computed range",
                  "of x values", prange(c(cl$from, cl$to))),
            call.=FALSE)
  #' integrate
  xx <- f$x
  yy <- f$y
  nn <- length(xx)
  Fx <- cumsum(yy * c(0, diff(xx)))
  #' normalise
  Fx <- Fx/Fx[nn]
  #' 
  FF <- approxfun(xx, Fx, method="linear", rule=2)
  return(FF)
}

