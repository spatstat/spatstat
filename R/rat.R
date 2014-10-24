#
#    rat.R
#
#   Ratio objects
#
#   Numerator and denominator are stored as attributes
#
#   $Revision: 1.5 $   $Date: 2014/10/24 00:22:30 $
#

rat <- function(ratio, numerator, denominator, check=TRUE) {
  if(check) {
    stopifnot(compatible(numerator, denominator))
    stopifnot(compatible(ratio, denominator))
  }
  attr(ratio, "numerator") <- numerator
  attr(ratio, "denominator") <- denominator
  class(ratio) <- c("rat", class(ratio))
  return(ratio)
}

print.rat <- function(x, ...) {
  NextMethod("print")
  cat("[Contains ratio information]\n")
  return(invisible(NULL))
}

compatible.rat <- function(A, B, ...) {
  NextMethod("compatible")
}

pool.rat <- function(...) {
  argh <- list(...)
  n <- narg <- length(argh)
  if(narg == 0) return(NULL)
  if(narg == 1) return(argh[[1]])
  #
  israt <- unlist(lapply(argh, inherits, what="rat"))
  if(any(bad <- !israt)) {
    nbad <- sum(bad)
    stop(paste(ngettext(nbad, "Argument", "Arguments"),
               commasep(which(bad)),
               ngettext(nbad, "does not", "do not"),
               "contain ratio (numerator/denominator) information"))
  }
  isfv <- unlist(lapply(argh, is.fv))
  if(!all(isfv))
    stop("All arguments must be fv objects")
  # extract
  template <- vanilla.fv(argh[[1]])
  Y <- lapply(argh, attr, which="numerator")
  X <- lapply(argh, attr, which="denominator")
  templateX <- vanilla.fv(X[[1]])
  templateY <- vanilla.fv(Y[[1]])
  # sum
  Add <- function(A,B){ force(A); force(B); eval.fv(A+B) }
  sumX <- Reduce(Add, X)
  sumY <- Reduce(Add, Y)
  attributes(sumX) <- attributes(templateX)
  attributes(sumY) <- attributes(templateY)
  # ratio-of-sums
  Ratio <- eval.fv(sumY/sumX)
  # variance calculation
  meanX <- eval.fv(sumX/n)
  meanY <- eval.fv(sumY/n)
  Square <- function(A) { force(A); eval.fv(A^2) }
  sumX2 <- Reduce(Add, lapply(X, Square))
  sumY2 <- Reduce(Add, lapply(Y, Square))
  varX   <- eval.fv((sumX2 - n * meanX^2)/(n-1))
  varY   <- eval.fv((sumY2 - n * meanY^2)/(n-1))
  Mul <- function(A,B){ force(A); force(B); eval.fv(A*B) }
  XY <- Map(Mul, X, Y)
  sumXY <- Reduce(Add, XY)
  covXY <- eval.fv((sumXY - n * meanX * meanY)/(n-1))
  # variance by delta method
  relvar <- eval.fv(pmax.int(0, varY/meanY^2 + varX/meanX^2
                            - 2 * covXY/(meanX * meanY)))
  Variance <- eval.fv(Ratio^2 * relvar/n)
  # two sigma CI
  hiCI <- eval.fv(Ratio + 2 * sqrt(Variance))
  loCI <- eval.fv(Ratio - 2 * sqrt(Variance))
  # relabel
  attributes(Ratio) <- attributes(Variance) <- attributes(template)
  Ratio <- prefixfv(Ratio,
                    tagprefix="pool",
                    descprefix="pooled ",
                    lablprefix="")
  Variance <- prefixfv(Variance,
                    tagprefix="var",
                    descprefix="delta-method variance estimate of ",
                    lablprefix="bold(var)~")
  attributes(hiCI) <- attributes(loCI) <-  attributes(template)
  hiCI <- prefixfv(hiCI,
                   tagprefix="hi",
                    descprefix="upper limit of two-sigma CI based on ",
                    lablprefix="bold(hi)~")
  loCI <- prefixfv(loCI,
                   tagprefix="lo",
                   descprefix="lower limit of two-sigma CI based on ",
                   lablprefix="bold(lo)~")
  #
  result <- Reduce(bind.fv, list(Ratio, Variance, hiCI, loCI))
  return(result)
}


