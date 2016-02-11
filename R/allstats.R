#
#
#   allstats.R
#
#   $Revision: 1.18 $   $Date: 2016/02/11 10:17:12 $
#
#
allstats <- function(pp, ..., dataname=NULL,verb=FALSE) {
#
# Function allstats --- to calculate the F, G, K, and J functions
# for an unmarked point pattern.
#
  verifyclass(pp,"ppp")
  if(is.marked(pp))
    stop("This function is applicable only to unmarked patterns.\n")

# estimate F, G and J 
  if(verb) cat("Calculating F, G, J ...")
  Jout <- do.call.matched(Jest,list(X=pp, ...))
  if(verb) cat("ok.\n")

# extract F, G and J
  Fout <- attr(Jout, "F")
  Gout <- attr(Jout, "G")
  attr(Jout, "F") <- NULL
  attr(Jout, "G") <- NULL
  fns <- list("F function"=Fout,
              "G function"=Gout,
              "J function"=Jout)

# compute second moment function K
  if(verb) cat("Calculating K function...")
  Kout <- do.call.matched(Kest, list(X=pp, ...))
  fns <- append(fns, list("K function"=Kout))
  if(verb) cat("done.\n")

# add title
  if(is.null(dataname))
    dataname <- short.deparse(substitute(pp))
  title <- paste("Four summary functions for ",
              	dataname,".",sep="")
  attr(fns, "title") <- title

#
  fns <- as.anylist(fns)
  return(fns)
}
