#
#    rat.R
#
#   Ratio objects
#
#   Numerator and denominator are stored as attributes
#
#   $Revision: 1.11 $   $Date: 2017/07/13 08:02:16 $
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

pool.rat <- local({

  Add <- function(A,B){ force(A); force(B); eval.fv(A+B, relabel=FALSE) }
  Square <- function(A) { force(A);         eval.fv(A^2, relabel=FALSE) }
  Mul <- function(A,B){ force(A); force(B); eval.fv(A*B, relabel=FALSE) }

  pool.rat <- function(..., weights=NULL, relabel=TRUE, variance=TRUE) {
    argh <- list(...)
    n <- narg <- length(argh)
    if(narg == 0) return(NULL)
    if(narg == 1) return(argh[[1]])
    ##
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
    ## extract
    template <- vanilla.fv(argh[[1]])
    Y <- lapply(argh, attr, which="numerator")
    X <- lapply(argh, attr, which="denominator")
    X <- do.call(harmonise, X)
    Y <- do.call(harmonise, Y)
    templateX <- vanilla.fv(X[[1]])
    templateY <- vanilla.fv(Y[[1]])
    ## compute products
    if(!is.null(weights)) {
      check.nvector(weights, narg, things="Functions")
      X <- Map(Mul, X, weights)
      Y <- Map(Mul, Y, weights)
    } 
    ## sum
    sumX <- Reduce(Add, X)
    sumY <- Reduce(Add, Y)
    attributes(sumX) <- attributes(templateX)
    attributes(sumY) <- attributes(templateY)
    ## ratio-of-sums
    Ratio <- eval.fv(sumY/sumX, relabel=FALSE)
    attributes(Ratio) <- attributes(template)
    ## variance calculation
    if(variance) { 
      meanX <- eval.fv(sumX/n, relabel=FALSE)
      meanY <- eval.fv(sumY/n, relabel=FALSE)
      sumX2 <- Reduce(Add, lapply(X, Square))
      sumY2 <- Reduce(Add, lapply(Y, Square))
      varX   <- eval.fv((sumX2 - n * meanX^2)/(n-1), relabel=FALSE)
      varY   <- eval.fv((sumY2 - n * meanY^2)/(n-1), relabel=FALSE)
      XY <- Map(Mul, X, Y)
      sumXY <- Reduce(Add, XY)
      covXY <- eval.fv((sumXY - n * meanX * meanY)/(n-1), relabel=FALSE)
      ## variance by delta method
      relvar <- eval.fv(pmax.int(0, varY/meanY^2 + varX/meanX^2
                                 - 2 * covXY/(meanX * meanY)),
		        relabel=FALSE)
      Variance <- eval.fv(Ratio^2 * relvar/n, relabel=FALSE)
      attributes(Variance) <- attributes(template)
      ## two sigma CI
      hiCI <- eval.fv(Ratio + 2 * sqrt(Variance), relabel=FALSE)
      loCI <- eval.fv(Ratio - 2 * sqrt(Variance), relabel=FALSE)
      attributes(hiCI) <- attributes(loCI) <-  attributes(template)
    }
    ## dress up
    if(relabel) {
      Ratio <- prefixfv(Ratio,
                        tagprefix="pool",
                        descprefix="pooled ",
                        lablprefix="")
      if(variance) {		      
        Variance <- prefixfv(Variance,
                             tagprefix="var",
                             descprefix="delta-method variance estimate of ",
                             lablprefix="bold(var)~")
        hiCI <- prefixfv(hiCI,
                         tagprefix="hi",
                         descprefix="upper limit of two-sigma CI based on ",
                         lablprefix="bold(hi)~")
        loCI <- prefixfv(loCI,
                         tagprefix="lo",
                         descprefix="lower limit of two-sigma CI based on ",
                         lablprefix="bold(lo)~")
      }
    }
    result <- if(!variance) Ratio else
              Reduce(bind.fv, list(Ratio, Variance, hiCI, loCI))
    return(result)
  }

  pool.rat
  
})

adjust.ratfv <- function(f, columns=fvnames(f, "*"), numfactor=1, denfactor=1) {
  stopifnot(is.fv(f))
  f[,columns] <- (numfactor/denfactor) * as.data.frame(f)[,columns]
  if(numfactor != 1 && !is.null(num <- attr(f, "numerator"))) {
    num[,columns] <- numfactor * as.data.frame(num)[,columns]
    attr(f, "numerator") <- num
  }	    
  if(denfactor != 1 && !is.null(den <- attr(f, "denominator"))) {
    den[,columns] <- denfactor * as.data.frame(den)[,columns]
    attr(f, "denominator") <- den
  }
  return(f)
}  

tweak.ratfv.entry <- function(x, ...) {
  # apply same tweak to function, numerator and denominator.
  x <- tweak.fv.entry(x, ...)
  if(!is.null(num <- attr(x, "numerator")))
    attr(x, "numerator") <- tweak.fv.entry(num, ...)
  if(!is.null(den <- attr(x, "denominator")))
    attr(x, "denominator") <- tweak.fv.entry(den, ...)
  return(x)
}

"[.rat" <- function(x, ...) {
   if(!is.fv(x)) stop("Not yet implemented for non-fv ratios")
   num <- attr(x, "numerator")
   den <- attr(x, "denominator")
   class(x) <- "fv"
   x <- x[...]
   den <- den[...]
   num <- num[...]
   attr(x, "numerator") <- num
   attr(x, "denominator") <- den
   class(x) <- c("rat", class(x))
   return(x)
}
  
