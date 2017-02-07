#
#
#   formulae.R
#
#   Functions for manipulating model formulae
#
#	$Revision: 1.25 $	$Date: 2017/02/07 07:22:47 $
#
#   identical.formulae()
#          Test whether two formulae are identical
#
#   termsinformula()
#          Extract the terms from a formula
#
#   sympoly()
#          Create a symbolic polynomial formula
#
#   polynom()
#          Analogue of poly() but without dynamic orthonormalisation
#
# -------------------------------------------------------------------
#	


identical.formulae <- function(x, y) {
  # workaround for bug in all.equal.formula in R 2.5.0
  if(is.null(y) && !is.null(x))
    return(FALSE)
  return(identical(all.equal(x,y), TRUE))
}

termsinformula <- function(x) {
  if(is.null(x)) return(character(0))
  if(class(x) != "formula")
    stop("argument is not a formula")
  attr(terms(x), "term.labels")
}

variablesinformula <- function(x) {
  if(is.null(x)) return(character(0))
  if(class(x) != "formula")
    stop("argument is not a formula")
  all.vars(as.expression(x))
}

offsetsinformula <- function(x) {
  if(is.null(x)) return(character(0))
  if(class(x) != "formula")
    stop("argument is not a formula")
  tums <- terms(x)
  offs <- attr(tums, "offset")
  if(length(offs) == 0) return(character(0))
  vars <- attr(tums, "variables")
  termnames <- unlist(lapply(vars, deparse))[-1L]
  termnames[offs]
}
  
lhs.of.formula <- function(x) {
  if(!inherits(x, "formula"))
    stop("x must be a formula")
  if(length(as.list(x)) == 3) {
    # formula has a response: return it
    return(x[[2L]])
  }
  return(NULL)
}

rhs.of.formula <- function(x, tilde=TRUE) {
  if(!inherits(x, "formula"))
    stop("x must be a formula")
  if(length(as.list(x)) == 3) {
    # formula has a response: strip it
    x <- x[-2L]
  }
  if(!tilde) # remove the "~"
    x <- x[[2L]]
  return(x)
}


sympoly <- function(x,y,n) {

   if(nargs()<2) stop("Degree must be supplied.")
   if(nargs()==2) n <- y
   eps <- abs(n%%1)
   if(eps > 0.000001 | n <= 0) stop("Degree must be a positive integer")
   
   x <- deparse(substitute(x))
   temp <- NULL
   left <- "I("
   rght <- ")"
   if(nargs()==2) {
	for(i in 1:n) {
		xhat <- if(i==1) "" else paste("^",i,sep="")
		temp <- c(temp,paste(left,x,xhat,rght,sep=""))
	}
   }
   else {
	y <- deparse(substitute(y))
	for(i in 1:n) {
		for(j in 0:i) {
			k <- i-j
			xhat <- if(k<=1) "" else paste("^",k,sep="")
			yhat <- if(j<=1) "" else paste("^",j,sep="")
			xbit <- if(k>0) x else ""
			ybit <- if(j>0) y else ""
			star <- if(j*k>0) "*" else ""
			term <- paste(left,xbit,xhat,star,ybit,yhat,rght,sep="")
			temp <- c(temp,term)
		}
	}
      }
   as.formula(paste("~",paste(temp,collapse="+")))
 }



expand.polynom <- local({
  powername <- function(x, n) {
    ifelse(n == 0, "",
           ifelse(n == 1,
                  x,
                  paste0(x, "^", n)))
  }
  power1name <- function(x, n) {
    px <- powername(x, n)
    ifelse(n <= 1, px, paste0("I", paren(px)))
  }
  power2name <- function(x, y, n, m) {
    ifelse(n == 0,
           power1name(y, m),
           ifelse(m == 0,
                  power1name(x, n),
                  paste0("I", paren(paste(powername(x, n),
                                          powername(y, m), sep="*")))))
  }

  haspolynom <- function(z) { 'polynom' %in% all.names(z) }

  fiddle <- function(f) {
    opname <- f[[1L]]
    if(identical(opname, as.name('I'))) {
      ## expressions enclosed in I() are protected
      return(f)
    }
    if(!identical(opname, as.name('polynom'))) {
      tbd <- unlist(lapply(f, haspolynom))
      if(any(tbd)) {
        ## descend recursively
        for(i in which(tbd)) 
          f[[i]] <- fiddle(f[[i]])
      }
      return(f)
    }
    ## polynom(..., d)
    n <- length(f)
    if(!(n %in% c(3,4)))
      stop("Syntax of polynom() call not understood")
    degree <- f[[n]]
    if (!is.numeric(degree) || length(degree) != 1 ||
        (degree%%1) != 0 || degree < 1) 
      stop("degree of polynomial should be a positive integer")
    if(n == 3) {
      ## polynom(x, d)
      xlang <- f[[2L]]
      xstring <- if(length(xlang) == 1) paste(xlang) else paren(format(xlang))
      xpowers <- power1name(xstring, 1:degree)
      xpolystring <- paste(xpowers, collapse=" + ")
      xpolylang <- as.formula(paste("~", xpolystring))[[2L]]
      return(xpolylang)
    } else if(n == 4) {
      ## polynom(x, y, d)
      xlang <- f[[2L]]
      ylang <- f[[3L]]
      xstring <- if(length(xlang) == 1) paste(xlang) else paren(format(xlang))
      ystring <- if(length(ylang) == 1) paste(ylang) else paren(format(ylang))
      mat <- matrix(, 1+degree, 1+degree)
      totdeg <- col(mat) - 1
      yd <- row(mat) - 1
      xd <- totdeg - yd
      xdeg <- xd[xd >= 0]
      ydeg <- yd[xd >= 0]
      xypowers <- power2name(xstring, ystring, xdeg, ydeg)[xdeg + ydeg > 0]
      xypolystring <- paste(xypowers, collapse=" + ")
      xypolylang <- as.formula(paste("~", xypolystring))[[2L]]
      return(xypolylang)
    }
  }

  expand.polynom <- function(f) {
    ## replaces polynom(...) by x + I(x^2) + ... inside a formula f
    g <- fiddle(f)
    environment(g) <- environment(f)
    return(g)
  }

  expand.polynom
})

can.be.formula <- function(x) {
  #' test whether x is a formula object
  if(inherits(x, "formula")) return(TRUE)
  #' or a character representation of a formula.
  if(!is.character(x)) return(FALSE)
  x <- paste(x, collapse=" ")
  if(length(grep("~", x)) == 0) return(FALSE)
  ok <- !inherits(try(as.formula(x), silent=TRUE), "try-error")
  return(ok)
}
