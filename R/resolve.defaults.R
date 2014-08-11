#
#   resolve.defaults.R
#
#  $Revision: 1.13 $ $Date: 2012/12/06 04:34:03 $
#
# Resolve conflicts between several sets of defaults
# Usage:
#     resolve.defaults(list1, list2, list3, .......)
# where the earlier lists have priority 
#
resolve.defaults <- function(..., .MatchNull=TRUE, .StripNull=FALSE) {
  # Each argument is a list. Append them.
  argue <- c(...)
  # is NULL a possible value?
  if(!.MatchNull) {
    isnul <- unlist(lapply(argue, is.null))
    argue <- argue[!isnul]
  }
  if(!is.null(nam <- names(argue))) {
    named <- nzchar(nam)
    arg.unnamed <- argue[!named]
    arg.named <-   argue[named]
    if(any(discard <- duplicated(names(arg.named)))) 
      arg.named <- arg.named[!discard]
    argue <- append(arg.unnamed, arg.named)
  }
  # should NULL become a missing argument?
  if(.StripNull) {
    isnull <- sapply(argue, is.null)
    argue <- argue[!isnull]
  }
  return(argue)
}

do.call.matched <- function(fun, arglist, funargs, extrargs=NULL, sieve=FALSE) {
  if(!is.function(fun) && !is.character(fun))
    stop("Internal error: wrong argument type in do.call.matched")
  if(is.character(fun)) {
    fname <- fun
    fun <- get(fname, mode="function")
    if(!is.function(fun))
      stop(paste("internal error: function", sQuote(fname), "not found",
                 sep=""))
  } 
  if(missing(funargs))
    funargs <- names(formals(fun))
  funargs <- c(funargs, extrargs)
  givenargs <- names(arglist)
  matched <- givenargs %in% funargs
  # apply 'fun' to matched arguments
  out <- do.call(fun, arglist[matched])
  # retain un-matched arguments?
  if(sieve)
    out <- list(result=out, otherargs=arglist[!matched])
  return(out)
}

resolve.1.default <- function(.A, ...) {
  res <- resolve.defaults(...)
  hit <- (names(res) == .A)
  if(!any(hit)) return(NULL)
  return(res[[min(which(hit))]])
}

# extract all the arguments that match '...' rather than a named argument

passthrough <- function(.Fun, ..., .Fname=NULL) {
  if(is.null(.Fname))
    .Fname <- deparse(substitute(.Fun))
  # make a fake call to the named function using the arguments provided
  cl <- eval(substitute(call(.Fname, ...)))
  # match the call to the function 
  mc <- match.call(.Fun, cl)
  # extract the arguments
  mcargs <- as.list(mc)[-1]
  # figure out which ones are actually formal arguments of the function
  nam <- names(formals(.Fun))
  nam <- setdiff(nam, "...")
  known <- names(mcargs) %in% nam
  # return the *other* arguments
  return(mcargs[!known])
}
