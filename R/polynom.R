#'
#'    polynom.R
#'
#'   $Revision: 1.1 $  $Date: 2017/01/02 09:48:36 $
#'

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
