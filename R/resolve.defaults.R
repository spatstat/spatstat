#
#   resolve.defaults.R
#
#  $Revision: 1.22 $ $Date: 2014/12/14 01:50:04 $
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

do.call.without <- function(fun, ..., avoid) {
  argh <- list(...)
  nama <- names(argh)
  if(!is.null(nama))
    argh <- argh[!(nama %in% avoid)]
  do.call(fun, argh)
}

do.call.matched <- function(fun, arglist, funargs,
                            extrargs=NULL,
                            matchfirst=FALSE,
                            sieve=FALSE,
                            skipargs=NULL) {
  if(!is.function(fun) && !is.character(fun))
    stop("Internal error: wrong argument type in do.call.matched")
  if(is.character(fun)) {
    fname <- fun
    fun <- get(fname, mode="function")
    if(!is.function(fun))
      stop(paste("internal error: function", sQuote(fname), "not found",
                 sep=""))
  }
  ## determine list of argument names to be matched
  if(missing(funargs))
    funargs <- names(formals(fun))
  funargs <- c(funargs, extrargs)
  funargs <- setdiff(funargs, skipargs)
  ## identify which arguments in the call actually match a formal argument
  givenargs <- names(arglist)
  matched <- givenargs %in% funargs
  # deem the first argument to be matched?
  if(matchfirst && !nzchar(givenargs[1]))
    matched[1] <- TRUE
  # apply 'fun' to matched arguments
  out <- do.call(fun, arglist[matched])
  # retain un-matched arguments?
  if(sieve)
    out <- list(result=out, otherargs=arglist[!matched])
  return(out)
}

## This function traps the colour arguments
## and converts to greyscale if required.

do.call.plotfun <- function(fun, arglist, ...) {
  if(spatstat.options("monochrome")) {
    keys <- names(arglist)
    if(!is.null(keys)) {
      cols <- nzchar(keys) & ((keys %in% c("border", "col", "fg", "bg")) |
                              (substr(keys, 1, 4) == "col."))
      if(any(cols))
        arglist[cols] <- lapply(arglist[cols], to.grey)
    }
  }
  do.call.matched(fun, arglist, ...)
}


resolve.1.default <- function(.A, ...) {
  if(is.character(.A)) {
    ## .A is the name of the parameter to be returned
    Aname <- .A
    res <- resolve.defaults(...)
  } else if(is.list(.A) && length(.A) == 1) {
    ## .A is a list giving the name and default value of the parameter
    Aname <- names(.A)
    res <- resolve.defaults(..., .A)
  } else stop("Unrecognised format for .A")
  hit <- (names(res) == Aname)
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

graphicsPars <- local({
  ## recognised additional arguments to image.default(), axis() etc
    PlotArgs <- c(
        "main", "asp", "sub", "axes", "ann",
        "cex", "font", 
        "cex.axis", "cex.lab", "cex.main", "cex.sub",
        "col.axis", "col.lab", "col.main", "col.sub",
        "font.axis", "font.lab", "font.main", "font.sub")
    
  TheTable <- 
    list(plot = PlotArgs,
         image = c(
           "main", "asp", "sub", "axes", "ann",
           "box", 
           "cex", "font", 
           "cex.axis", "cex.lab", "cex.main", "cex.sub",
           "col.axis", "col.lab", "col.main", "col.sub",
           "font.axis", "font.lab", "font.main", "font.sub"),
         axis = c(
           "cex", 
           "cex.axis", "cex.lab",
           "col.axis", "col.lab",
           "font.axis", "font.lab",
           "mgp", "xaxp", "yaxp", "tck", "tcl", "las", "fg", "xpd"),
         owin = c(
           "sub",
           "cex", "font", "col",
           "border", "box", 
           "cex.main", "cex.sub",
           "col.main", "col.sub",
           "font.main", "font.sub",
           "xaxs", "yaxs"),
         lines = c("lwd", "lty", "col", "lend", "ljoin", "lmitre"),
         symbols = c(PlotArgs, "fg", "bg")
         )

    TheTable$ppp <- unique(c(TheTable$owin,
                             TheTable$symbols,
                             "pch", "cex", "lty", "lwd",
                             "etch"))

  graphicsPars <- function(key) {
    n <- pmatch(key, names(TheTable))
    if(is.na(n)) return(NULL)
    return(TheTable[[n]])
  }

  graphicsPars
})
