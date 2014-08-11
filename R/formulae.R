#
#
#   formulae.S
#
#   Functions for manipulating model formulae
#
#	$Revision: 1.20 $	$Date: 2014/04/14 08:05:01 $
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

# new generic

"formula<-" <- function(x, ..., value) {
  UseMethod("formula<-")
}


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
  termnames <- unlist(lapply(vars, deparse))[-1]
  termnames[offs]
}
  
lhs.of.formula <- function(x) {
  if(!inherits(x, "formula"))
    stop("x must be a formula")
  if(length(as.list(x)) == 3) {
    # formula has a response: return it
    return(x[[2]])
  }
  return(NULL)
}

rhs.of.formula <- function(x, tilde=TRUE) {
  if(!inherits(x, "formula"))
    stop("x must be a formula")
  if(length(as.list(x)) == 3) {
    # formula has a response: strip it
    x <- x[-2]
  }
  if(!tilde) # remove the "~"
    x <- x[[2]]
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


polynom <- function(x, ...) {
  rest <- list(...)
  # degree not given
  if(length(rest) == 0)
    stop("degree of polynomial must be given")
  #call with single variable and degree
  if(length(rest) == 1) {
    degree <- ..1
    if((degree %% 1) != 0 || length(degree) != 1 || degree < 1)
      stop("degree of polynomial should be a positive integer")

    # compute values
    result <- outer(x, 1:degree, "^")

    # compute column names - the hard part !
    namex <- deparse(substitute(x))
    # check whether it needs to be parenthesised
    if(!is.name(substitute(x))) 
      namex <- paste("(", namex, ")", sep="")
    # column names
    namepowers <- if(degree == 1) namex else 
                       c(namex, paste(namex, "^", 2:degree, sep=""))
    namepowers <- paste("[", namepowers, "]", sep="")
    # stick them on
    dimnames(result) <- list(NULL, namepowers)
    return(result)
  }
  # call with two variables and degree
  if(length(rest) == 2) {

    y <- ..1
    degree <- ..2

    # list of exponents of x and y, in nice order
    xexp <- yexp <- numeric()
    for(i in 1:degree) {
      xexp <- c(xexp, i:0)
      yexp <- c(yexp, 0:i)
    }
    nterms <- length(xexp)
    
    # compute 

    result <- matrix(, nrow=length(x), ncol=nterms)
    for(i in 1:nterms) 
      result[, i] <- x^xexp[i] * y^yexp[i]

    #  names of these terms
    
    namex <- deparse(substitute(x))
    # namey <- deparse(substitute(..1)) ### seems not to work in R
    zzz <- as.list(match.call())
    namey <- deparse(zzz[[3]])

    # check whether they need to be parenthesised
    # if so, add parentheses
    if(!is.name(substitute(x))) 
      namex <- paste("(", namex, ")", sep="")
    if(!is.name(zzz[[3]])) 
      namey <- paste("(", namey, ")", sep="")

    nameXexp <- c("", namex, paste(namex, "^", 2:degree, sep=""))
    nameYexp <- c("", namey, paste(namey, "^", 2:degree, sep=""))

    # make the term names
       
    termnames <- paste(nameXexp[xexp + 1],
                       ifelse(xexp > 0 & yexp > 0, ".", ""),
                       nameYexp[yexp + 1],
                       sep="")
    termnames <- paste("[", termnames, "]", sep="")

    dimnames(result) <- list(NULL, termnames)
    # 
    return(result)
  }
  stop("Can't deal with more than 2 variables yet")
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

  fiddle <- function(f) {
    opname <- f[[1]]
    if(identical(opname, as.name('I'))) {
      ## expressions enclosed in I() are protected
      return(f)
    }
    if(!identical(opname, as.name('polynom'))) {
      tbd <- unlist(lapply(f, function(z) { 'polynom' %in% all.names(z) }))
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
      xlang <- f[[2]]
      xstring <- if(length(xlang) == 1) paste(xlang) else paren(format(xlang))
      xpowers <- power1name(xstring, 1:degree)
      xpolystring <- paste(xpowers, collapse=" + ")
      xpolylang <- as.formula(paste("~", xpolystring))[[2]]
      return(xpolylang)
    } else if(n == 4) {
      ## polynom(x, y, d)
      xlang <- f[[2]]
      ylang <- f[[3]]
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
      xypolylang <- as.formula(paste("~", xypolystring))[[2]]
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
