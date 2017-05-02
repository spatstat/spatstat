#
#      envelopeArray.R
#
#   $Revision$   $Date$
#
#
                                  
envelopeArray <- function(X, fun, ...,
                          dataname=NULL,verb=FALSE,reuse=TRUE) {
#'
  if(is.null(dataname))
    dataname <- short.deparse(substitute(X))

#' determine function name
  f.is.name <- is.name(substitute(fun))
  fname <-
    if(f.is.name)
      paste(as.name(substitute(fun)))
    else if(is.character(fun))
      fun
    else sQuote("fun") 

#' determine function to be called

  if(is.character(fun)) {
    fun <- get(fun, mode="function")
  } else if(!is.function(fun)) 
    stop(paste(sQuote("fun"), "should be a function or a character string"))

#' Apply function to data pattern, to test it
#' and to determine array dimensions, margin labels etc.

  fX <- do.call.matched(fun, append(list(X), list(...)), matchfirst=TRUE)
  if(!inherits(fX, "fasp"))
     stop("function did not return an object of class 'fasp'")

  d <- dim(fX)
  witch <- matrix(1:prod(d), nrow=d[1L], ncol=d[2L],
                  dimnames=dimnames(fX))

#' make function that extracts [i,j] entry of result

   ijfun <- function(X, ..., i=1, j=1, expectdim=d) {
     fX <- fun(X, ...)
     if(!inherits(fX, "fasp"))
       stop("function did not return an object of class 'fasp'")
     if(!all(dim(fX) == expectdim))
       stop("function returned an array with different dimensions")
     return(fX[i,j])
   }
   
  # ------------ start computing -------------------------------  
  if(reuse) {
    L <- do.call(spatstat::envelope,
                 resolve.defaults(
                                  list(X, fun=ijfun),
                                  list(internal=list(eject="patterns")),
                                  list(...),
                                  list(verbose=verb)))
    intern <- attr(L, "internal")
  } else intern <- L <- NULL

  # compute function array and build up 'fasp' object
  fns  <- list()
  k   <- 0

  for(i in 1:nrow(witch)) {
    for(j in 1:ncol(witch)) {
      if(verb) cat("i =",i,"j =",j,"\n")
      currentfv <- 
        do.call(spatstat::envelope,
                resolve.defaults(
                                 list(X, ijfun),
                                 list(simulate=L, internal=intern),
                                 list(verbose=FALSE),
                                 list(...),
                                 list(Yname=dataname),
				 list(i=i, j=j)))
      k <- k+1
      fns[[k]] <- as.fv(currentfv)
    }
  }

  # wrap up into 'fasp' object
  title <- paste("array of envelopes of", fname,
                 "for", dataname)
  
  rslt <- fasp(fns, which=witch,
               formulae=NULL,
               dataname=dataname,
               title=title,
               checkfv=FALSE)
  return(rslt)
}

