#
#    rat.R
#
#   Ratio objects
#
#   Numerator and denominator are stored as attributes
#
#   $Revision: 1.13 $   $Date: 2020/11/30 09:43:44 $
#

rat <- function(ratio, numerator, denominator, check=TRUE) {
  if(check) {
    stopifnot(compatible(numerator, denominator))
    stopifnot(compatible(ratio, denominator))
  }
  attr(ratio, "numerator") <- numerator
  attr(ratio, "denominator") <- denominator
  class(ratio) <- unique(c("rat", class(ratio)))
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
   class(x) <- unique(c("rat", class(x)))
   return(x)
}
  
