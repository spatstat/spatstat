#
#
#	classes.S
#
#	$Revision: 1.7 $	$Date: 2006/10/09 03:38:14 $
#
#	Generic utilities for classes
#
#
#--------------------------------------------------------------------------

verifyclass <- function(X, C, N=deparse(substitute(X)), fatal=TRUE) {
  if(!inherits(X, C)) {
    if(fatal) {
        gripe <- paste("argument", sQuote(N),
                       "is not of class", sQuote(C))
	stop(gripe)
    } else 
	return(FALSE)
  }
  return(TRUE)
}

#--------------------------------------------------------------------------

checkfields <- function(X, L) {
	  # X is a list, L is a vector of strings
	  # Checks for presence of field named L[i] for all i
	return(all(!is.na(match(L,names(X)))))
}

getfields <- function(X, L, fatal=TRUE) {
	  # X is a list, L is a vector of strings
	  # Extracts all fields with names L[i] from list X
	  # Checks for presence of all desired fields
	  # Returns the sublist of X with fields named L[i]
	absent <- is.na(match(L, names(X)))
	if(any(absent)) {
		gripe <- paste("Needed the following components:",
				paste(L, collapse=", "),
				"\nThese ones were missing: ",
				paste(L[absent], collapse=", "))
		if(fatal)
			stop(gripe)
		else 
			warning(gripe)
	} 
	return(X[L[!absent]])
}



