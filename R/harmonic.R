#
#
#   harmonic.R
#
#	$Revision: 1.2 $	$Date: 2004/01/07 08:57:39 $
#
#   harmonic()
#          Analogue of polynom() for harmonic functions only
#
# -------------------------------------------------------------------
#	

harmonic <- function(x,y,n) {
  if(missing(n))
    stop("the order n must be specified")
  n <- as.integer(n)
  if(is.na(n) || n <= 0)
    stop("n must be a positive integer")

  if(n > 3)
    stop("Sorry, harmonic() is not implemented for degree > 3")

  namex <- deparse(substitute(x))
  namey <- deparse(substitute(y))
  if(!is.name(substitute(x))) 
      namex <- paste("(", namex, ")", sep="")
  if(!is.name(substitute(y))) 
      namey <- paste("(", namey, ")", sep="")
  
  switch(n,
         {
           result <- cbind(x, y)
           names <- c(namex, namey)
         },
         {
           result <- cbind(x, y,
                           x*y, x^2-y^2)
           names <- c(namex, namey,
                      paste("(", namex, ".", namey, ")", sep=""),
                      paste("(", namex, "^2-", namey, "^2)", sep=""))
         },
         {
           result <- cbind(x, y,
                           x * y, x^2-y^2, 
                           x^3 - 3 * x * y^2, y^3 - 3 * x^2 * y)
           names <- c(namex, namey,
                      paste("(", namex, ".", namey, ")", sep=""),
                      paste("(", namex, "^2-", namey, "^2)", sep=""),
                      paste("(", namex, "^3-3", namex, ".", namey, "^2)",
                            sep=""),
                      paste("(", namey, "^3-3", namex, "^2.", namey, ")",
                            sep="")
                      )
         }
         )
  dimnames(result) <- list(NULL, names)
  return(result)
}
