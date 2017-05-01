#
#      alltypes.R
#
#   $Revision: 1.34 $   $Date: 2016/10/23 10:36:58 $
#
#
                                  
alltypes <- function(X, fun="K", ...,
                     dataname=NULL,verb=FALSE,envelope=FALSE,reuse=TRUE) {
#
# Function 'alltypes' --- calculates a summary function for
# each type, or each pair of types, in a multitype point pattern
#
  if(is.ppp(X)) classname <- "ppp" else
  if(is.lpp(X)) classname <- "lpp" else
  stop("X should be a ppp or lpp object")
  
  if(is.null(dataname))
    dataname <- short.deparse(substitute(X))

# --------------------------------------------------------------------  
# First inspect marks

  if(!is.marked(X)) {
    nmarks <- 0
    marklabels <- ""
  } else {
    if(!is.multitype(X))
      stop("the marks must be a factor")
    # ensure type names are parseable (for mathematical labels)
    levels(marks(X)) <- make.parseable(levels(marks(X)))
    mks <- marks(X)
    ma <- levels(mks)
    nmarks <- length(ma)
    marklabels <- paste(ma)
  }

# ---------------------------------------------------------------------
# determine function name

  f.is.name <- is.name(substitute(fun))
  fname <-
    if(f.is.name)
      paste(as.name(substitute(fun)))
    else if(is.character(fun))
      fun
    else sQuote("fun") 

# ---------------------------------------------------------------------
# determine function to be called
  
  if(is.function(fun)) {
    estimator <- fun
  } else if(is.character(fun)) {
    # First try matching one of the standard abbreviations K, G etc
    estimator <- getSumFun(fun, classname, (nmarks > 0), fatal=FALSE)
    if(is.null(estimator))
      estimator <- get(fun, mode="function")
  } else 
      stop(paste(sQuote("fun"), "should be a function or a character string"))
  
# ------------------------------------------------------------------  
# determine how the function shall be called.
#
  indices.expected <- sum(c("i", "j") %in% names(formals(estimator)))

  apply.to.split   <- (indices.expected == 0 && nmarks > 1)
  if(apply.to.split)
    ppsplit <- split(X)
  
# --------------------------------------------------------------------  
# determine array dimensions and margin labels
  witch <-
    if(nmarks == 0)
      matrix(1L, nrow=1L, ncol=1L, dimnames=list("",""))
    else if (nmarks == 1) 
      matrix(1L, nrow=1L, ncol=1L, dimnames=list(marklabels, marklabels))
    else if(indices.expected != 2)
      matrix(1L:nmarks, nrow=nmarks, ncol=1L,
             dimnames=list(marklabels, ""))
    else 
      matrix(1L:(nmarks^2),ncol=nmarks,nrow=nmarks, byrow=TRUE,
             dimnames=list(marklabels, marklabels))

  # ------------ start computing -------------------------------  
  # if computing envelopes, first generate simulated patterns
  # using undocumented feature of envelope()
  if(envelope && reuse) {
    L <- do.call(spatstat::envelope,
                 resolve.defaults(
                                  list(X, fun=estimator),
                                  list(internal=list(eject="patterns")),
                                  list(...),
				  switch(1L+indices.expected,
                                          NULL,
                                          list(i=ma[1L]),
                                          list(i=ma[1L], j=ma[2L]),
                                          NULL),
                                  list(verbose=verb)))
    intern <- attr(L, "internal")
  } else intern <- L <- NULL

  # compute function array and build up 'fasp' object
  fns  <- list()
  k   <- 0

  for(i in 1L:nrow(witch)) {
    Y <- if(apply.to.split) ppsplit[[i]] else X
    for(j in 1L:ncol(witch)) {
      if(verb) cat("i =",i,"j =",j,"\n")
      currentfv <- 
        if(!envelope) 
          switch(1L+indices.expected,
                 estimator(Y, ...),
                 estimator(Y, i=ma[i], ...),
                 estimator(Y, i=ma[i], j=ma[j], ...))
        else
          do.call(spatstat::envelope,
                  resolve.defaults(
                                   list(Y, estimator),
                                   list(simulate=L, internal=intern),
                                   list(verbose=FALSE),
                                   list(...),
                                   list(Yname=dataname),
                                   switch(1L+indices.expected,
                                          NULL,
                                          list(i=ma[i]),
                                          list(i=ma[i], j=ma[j]),
                                          NULL)))
      k <- k+1
      fns[[k]] <- as.fv(currentfv)
    }
  }

  # wrap up into 'fasp' object
  title <- paste(if(nmarks > 1) "array of " else NULL,
                 if(envelope) "envelopes of " else NULL,
                 fname,
                 if(nmarks <= 1) " function " else " functions ",
                 "for ", dataname, ".", sep="")
  
  rslt <- fasp(fns, which=witch,
               formulae=NULL,
               dataname=dataname,
               title=title,
               checkfv=FALSE)
  return(rslt)
}

# Lookup table for standard abbreviations of functions

getSumFun <- local({

  ftable <-
  rbind(
        data.frame(class="ppp", marked=FALSE,
                   abbrev=c("F", "G", "J", "K", "L", "pcf"),
                   full=c("Fest", "Gest", "Jest", "Kest", "Lest", "pcf"),
                   stringsAsFactors=FALSE),
        data.frame(class="ppp", marked=TRUE,
                   abbrev=c("F", "G", "J", "K", "L", "pcf"),
                   full=  c("Fest",
                     "Gcross", "Jcross", "Kcross", "Lcross",
                     "pcfcross"),
                   stringsAsFactors=FALSE),
        data.frame(class="lpp", marked=FALSE,
                   abbrev=c("K", "pcf"),
                   full=c("linearK", "linearpcf"),
                   stringsAsFactors=FALSE),
        data.frame(class="lpp", marked=TRUE,
                   abbrev=c("K", "pcf"),
                   full=c("linearKcross", "linearpcfcross"),
                   stringsAsFactors=FALSE)
        )

  getfun <- function(abbreviation, classname, ismarked, fatal=TRUE) {
    matches <- with(ftable,
                    which(abbrev == abbreviation &
                          class == classname &
                          marked == ismarked))
    if(length(matches) == 0) {
      if(!fatal)
        return(NULL)
      stop(paste("No match to function abbreviation",
                 sQuote(abbreviation),
                 "for class",
                 sQuote(classname)))
    }
    if(length(matches) > 1)
      stop("Ambiguous function name")
    fullname <- ftable$full[matches]
    get(fullname, mode="function")
  }

  getfun
})


