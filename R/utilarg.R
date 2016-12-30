#'
#'   utilarg.R
#'
#'   Utilities for checking/handling arguments
#'
#'  $Revision: 1.2 $  $Date: 2016/12/30 03:32:21 $
#'

"%orifnull%" <- function(a, b) {
  if(!is.null(a)) return(a)
  # b is evaluated only now
  return(b)
}

check.nvector <- function(v, npoints=NULL, fatal=TRUE, things="data points",
                          naok=FALSE, warn=FALSE, vname, oneok=FALSE) {
  # vector of numeric values for each point/thing
  if(missing(vname))
    vname <- sQuote(deparse(substitute(v)))
  whinge <- NULL
  nv <- length(v)
  if(!is.numeric(v))
    whinge <- paste(vname, "is not numeric")
  else if(!is.atomic(v) || !is.null(dim(v)))  # vector with attributes
    whinge <- paste(vname, "is not a vector")
  else if(!(is.null(npoints) || (nv == npoints)) &&
          !(oneok && nv == 1)) 
    whinge <- paste("The length of", vname,
                    paren(paste0("=", nv)), 
                    "should equal the number of", things,
                    paren(paste0("=", npoints)))
  else if(!naok && anyNA(v))
    whinge <- paste("Some values of", vname, "are NA or NaN")
  #
  if(!is.null(whinge)) {
    if(fatal) stop(whinge)
    if(warn) warning(whinge)
    ans <- FALSE
    attr(ans, "whinge") <- whinge
    return(ans)
  }
  return(TRUE)
}

check.nmatrix <- function(m, npoints=NULL, fatal=TRUE, things="data points",
                          naok=FALSE, squarematrix=TRUE, matchto="nrow",
                          warn=FALSE) {
  ## matrix of values for each thing or each pair of things
  mname <- sQuote(deparse(substitute(m)))
  whinge <- NULL
  if(!is.matrix(m))
    whinge <- paste(mname, "should be a matrix")
  else if(squarematrix && (nrow(m) != ncol(m)))
    whinge <- paste(mname, "should be a square matrix")
  else if(!naok && anyNA(m))
    whinge <- paste("Some values of", mname, "are NA or NaN")
  else if(!is.null(npoints)) {
    if(matchto=="nrow" && nrow(m) != npoints)
      whinge <- paste("Number of rows in", mname,
                      paren(paste0("=", nrow(m))),
                      "does not match number of", things,
                      paren(paste0("=", npoints)))
    else if(matchto=="ncol" && ncol(m) != npoints)
      whinge <- paste("Number of columns in", mname,
                      paren(paste0("=", ncol(m))),
                      "does not match number of", things,
                      paren(paste0("=", npoints)))
  }
  ##
  if(!is.null(whinge)) {
    if(fatal) stop(whinge)
    if(warn) warning(whinge)
    return(FALSE)
  }
  return(TRUE)
}

check.named.vector <- function(x, nam, context="", namopt=character(0),
                               onError=c("fatal", "null")) {
  xtitle <- deparse(substitute(x))
  onError <- match.arg(onError)
  problem <- check.named.thing(x, nam, namopt, sQuote(xtitle),
                               is.numeric(x), "vector", context,
                               fatal=(onError == "fatal"))
  if(length(problem) > 0 && onError == "null")
    return(NULL)
  opt <- namopt %in% names(x)
  return(x[c(nam, namopt[opt])])
}

check.named.list <- function(x, nam, context="", namopt=character(0),
                               onError=c("fatal", "null")) {
  xtitle <- deparse(substitute(x))
  onError <- match.arg(onError)
  problem <- check.named.thing(x, nam, namopt, sQuote(xtitle),
                               is.list(x), "list", context,
                               fatal=(onError == "fatal"))
  if(length(problem) > 0 && onError == "null")
    return(NULL)
  opt <- namopt %in% names(x)
  return(x[c(nam, namopt[opt])])
}

check.named.thing <- function(x, nam, namopt=character(0), xtitle=NULL,
                              valid=TRUE, type="object", context="",
                              fatal=TRUE) {
  if(is.null(xtitle))
    xtitle <- sQuote(deparse(substitute(x)))
  # check whether names(x) contains all obligatory names 'nam'
  # and possibly some of the optional names 'namopt'
  namesx <- names(x)
  omitted <- !(nam %in% namesx)
  foreign <- !(namesx %in% c(nam, namopt))
  if(valid && !any(omitted) && !any(foreign))
    return(character(0))
  # some condition violated
  if(nzchar(context))
    xtitle <- paste(context, xtitle)
  whinge <- paste(xtitle,
                  "must be a named", paste(type, ",", sep=""),
                  "with components", commasep(nam))
  if(length(namopt) > 0)
    whinge <- paste(whinge, paren(paste("and optionally", commasep(namopt))))
  if(any(omitted)) {
    grizzle <- paste(ngettext(sum(omitted), "parameter", "parameters"),
                     commasep(nam[omitted]),
                     "omitted")
    whinge <- paste(whinge, grizzle, sep="; ")
  }
  if(any(foreign)) {
    grizzle <- paste(ngettext(sum(foreign), "component", "components"),
                     commasep(namesx[foreign]),
                     "not recognised")
    whinge <- paste(whinge, grizzle, sep="; ")
  }
  if(fatal)
    stop(whinge, call.=FALSE)
  return(whinge)
}


validposint <- function(n, caller, fatal=TRUE) {
  nname <- deparse(substitute(n))
  if(length(n) != 1 || n != round(n) || n <=0) {
    if(!fatal)
      return(FALSE)
    prefix <- if(!missing(caller)) paste("In ", caller, ", ", sep="") else NULL
    stop(paste0(prefix, nname, "should be a single positive integer"),
         call.=FALSE)
  }
  return(TRUE)
}


# errors and checks

forbidNA <- function(x, context="", xname, fatal=TRUE, usergiven=TRUE) {
  if(missing(xname)) xname <- sQuote(deparse(substitute(x)))
  if(anyNA(x)) {
    if(usergiven) {
      # argument came from user
      offence <- ngettext(length(x), "be NA", "contain NA values")
      whinge <- paste(context, xname, "must not", offence)
    } else {
      # argument was computed internally
      violates <- ngettext(length(x), "is NA", "contains NA values")
      whinge <- paste(context, xname, violates)
    }
    if(fatal) stop(whinge, call.=FALSE)
    warning(whinge, call.=FALSE)
    return(FALSE)
  }
  return(TRUE)
}

check.finite <- function(x, context="", xname, fatal=TRUE, usergiven=TRUE) {
  if(missing(xname)) xname <- sQuote(deparse(substitute(x)))
  forbidNA(x, context, xname, fatal=fatal, usergiven=usergiven)
  if(any(!is.finite(x))) {
    if(usergiven) {
      # argument came from user
      oblige <- ngettext(length(x),
                         "be a finite value", "contain finite values")
      whinge <- paste(context, xname, "must", oblige)
    } else {
      # argument was computed internally
      violates <- ngettext(length(x),
                           "is not finite", "contains non-finite values")
      whinge <- paste(context, xname, violates)
    }
    if(fatal) stop(whinge, call.=FALSE)
    warning(whinge, call.=FALSE)
    return(FALSE)
  }
  return(TRUE)
}

complaining <- function(whinge, fatal=FALSE, value=NULL) {
  if(fatal) stop(whinge, call.=FALSE)
  warning(whinge, call.=FALSE)
  return(value)
}

check.1.real <- function(x, context="", fatal=TRUE) {
  xname <- deparse(substitute(x))
  if(!is.numeric(x) || length(x) != 1) {
    whinge <-  paste(sQuote(xname), "should be a single number")
    if(nzchar(context)) whinge <- paste(context, whinge)
    return(complaining(whinge, fatal=fatal, value=FALSE))
  }
  return(TRUE)
}

check.1.integer <- function(x, context="", fatal=TRUE) {
  xname <- deparse(substitute(x))
  if(!is.numeric(x) || length(x) != 1 || !is.finite(x) || x %% 1 != 0) {
    whinge <-  paste(sQuote(xname), "should be a single finite integer")
    if(nzchar(context)) whinge <- paste(context, whinge)
    return(complaining(whinge, fatal=fatal, value=FALSE))
  }
  return(TRUE)
}

check.1.string <- function(x, context="", fatal=TRUE) {
  xname <- deparse(substitute(x))
  if(!is.character(x) || length(x) != 1) {
    whinge <-  paste(sQuote(xname), "should be a single character string")
    if(nzchar(context)) whinge <- paste(context, whinge)
    return(complaining(whinge, fatal=fatal, value=FALSE))
  }
  return(TRUE)
}

explain.ifnot <- function(expr, context="") {
  ex <- deparse(substitute(expr))
  ans <- expr
  if(!(is.logical(ans) && length(ans) == 1 && ans))
    stop(paste(context, "it must be TRUE that", sQuote(ex)), call.=FALSE)
}

warn.ignored.args <- function(..., context=NULL) {
  if((narg <- length(list(...))) > 0) {
    whinge <- paste(narg, "unrecognised",
                    ngettext(narg, "argument was", "arguments were"),
                    "ignored")
    if(!is.null(context)) whinge <- paste(context, whinge)
    warning(context)
  }
}

trap.extra.arguments <- function(..., .Context="", .Fatal=FALSE) {
  z <- list(...)
  if((narg <- length(z)) == 0) return(FALSE)
  nama <- names(z)
  named <- nzchar(nama)
  whinge <- paste(.Context, ":", sep="")
  if(any(named)) {
    # some arguments are named: ensure all are named
    nama <- sQuote(nama)
    if(!all(named)) 
      nama[!named] <- paste("[Arg", 1:length(nama), ,"]", sep="")[!named]
    whinge <- paste(whinge,
                    "unrecognised",
                    ngettext(narg, "argument", "arguments"),
                    commasep(nama),
                    ngettext(narg, "was", "were"), "ignored")
  } else {
    # all arguments unnamed
    whinge <- paste(whinge, 
                    narg, "unrecognised",
                    ngettext(narg, "argument was", "arguments were"),
                    "ignored")   
  }
  if(.Fatal) stop(whinge, call.=FALSE) else warning(whinge, call.=FALSE)
  return(TRUE)
}


## replace recognised keywords by other keywords
mapstrings <- function(x, map=NULL) {
  if(is.null(map)) return(x)
  x <- as.character(x)
  from <- names(map)
  map <- as.character(map)
  if(sum(nzchar(from)) != length(map))
    stop("input names are missing in map", call.=FALSE)
  if(any(duplicated(from)))
    stop("input names are duplicated in map", call.=FALSE)
  i <- match(x, from)
  hit <- !is.na(i)
  x[hit] <- map[i[hit]]
  return(x)
}

there.is.no.try <- function(...) {
  #' do, or do not
  y <- try(..., silent=TRUE)
  if(inherits(y, "try-error")) return(NULL)
  return(y)
}

dont.complain.about <- function(...) {
  #' prevents code checkers complaining about 'unused variables'
  #' Typically needed where the variables in question
  #' are referenced in an expression that will be evaluated elsewhere. 
  return(invisible(NULL))
}


matchNameOrPosition <- function(expected, avail) {
  if(length(avail) < length(expected))
    stop("Not enough arguments to match", call.=FALSE)
  j <- match(expected, avail)
  if(!anyNA(j)) return(j)
  everything <- seq_along(avail)
  for(k in seq_along(expected)) {
    if(is.na(j[k]))
      j[k] <- min(setdiff(everything, j[-k]))
  }
  return(j)
}

badprobability <- function(x, NAvalue=NA) {
  ifelse(is.na(x), NAvalue, !is.finite(x) | x < 0 | x > 1)
}

# test for equivalence of two functions 
samefunction <- function(f, g) {
  identical(deparse(f), deparse(g))
}

#' .................. calls and expressions ..................


fakecallstring <- function(fname, parlist) {
  cl <- do.call(call, append(list(name = fname), parlist))
  return(format(cl))
}

dotexpr.to.call <- function(expr, dot="funX", evaluator="eval.fv") {
  # convert an expression into a function call
  # replacing "." by the specified variable 
  stopifnot(is.expression(expr))
  aa <- substitute(substitute(ee, list(.=as.name(d))),
                   list(ee=expr, d=dot))
  bb <- eval(parse(text=deparse(aa)))
  cc <- as.call(bb)
  cc[[1]] <- as.name(evaluator)
  return(cc)
}

inject.expr <- function(base, expr) {
  ## insert an expression inside a call and parse it
  txt <- sub(".", as.character(expr), as.character(base), fixed=TRUE)
  parse(text=txt)
}

  
## Match variable names to objects in 'data' list or environment
getdataobjects <- function(nama, envir, datalist=NULL, fatal=FALSE) {
  if(is.null(nama)) return(NULL)
  stopifnot(is.character(nama))
  n <- length(nama)
  y <- vector(mode="list", length=n)
  names(y) <- nama
  if(!is.null(datalist)) {
    hit <- nama %in% names(datalist)
    if(any(hit))
      y[hit] <- as.list(datalist)[nama[hit]]
    external <- unlist(lapply(y, is.null))
  } else external <- rep(TRUE, n)
  y[external] <- mget(nama[external], envir=envir,
                    ifnotfound=list(NULL), inherits=TRUE)
  if(fatal && any(bad <- unlist(lapply(y, is.null)))) {
    nbad <- sum(bad)
    stop(paste(ngettext(nbad, "Covariate", "Covariates"),
               commasep(sQuote(nama[bad])),
               ngettext(nbad, "was not found", "were not found")),
         call.=FALSE)
  }
  names(y) <- nama
  attr(y, "external") <- external
  return(y)
}

