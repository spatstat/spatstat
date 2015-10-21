#
# interactions.R
#
# Works out which interaction is in force for a given point pattern
#
#  $Revision: 1.16 $  $Date: 2015/10/21 09:06:57 $
#
#
impliedpresence <- function(tags, formula, df, extranames=character(0)) {
  # Determines, for each row of the data frame df,
  # whether the variable called tags[j] is required in the formula
  stopifnot(is.data.frame(df))
  stopifnot(inherits(formula, "formula"))
  stopifnot(is.character(tags))
  stopifnot(is.character(extranames))
#  allvars <- variablesinformula(formula)
  if(any(tags %in% names(df)))
    stop(paste(sQuote("tags"),
               "conflicts with the name of a column of",
               sQuote("df")))
  if(any(extranames %in% names(df)))
    stop(paste(sQuote("extranames"),
               "conflicts with the name of a column of",
               sQuote("df")))
  # answer is a matrix 
  nvars <- length(tags)
  nrows <- nrow(df)
  answer <- matrix(TRUE, nrows, nvars)
  # expand data frame with zeroes for each tags and extranames
  for(v in unique(c(tags, extranames)))
    df[ , v] <- 0
  # loop
  for(i in seq(nrow(df))) {
    # make a fake data frame for the formula
    # using the data frame entries from row i
    # (includes 0 values for all other variables)
    pseudat <- df[i, , drop=FALSE]
    # use this to construct a fake model matrix
    mof0 <- model.frame(formula, pseudat)
    mom0 <- model.matrix(formula, mof0)
    for(j in seq(nvars)) {
      # Reset the variable called tags[j] to 1
      pseudatj <- pseudat
      pseudatj[ , tags[j]] <- 1
      # Now create the fake model matrix
      mofj <- model.frame(formula, pseudatj)
      momj <- model.matrix(formula, mofj)
      # Compare the two matrices
      answer[i,j] <- any(momj != mom0)
    }
  }
  return(answer)
}

active.interactions <- function(object) {
  stopifnot(inherits(object, "mppm"))
  interaction <- object$Inter$interaction
  iformula <- object$iformula
  nenv <- new.env()
  environment(iformula) <- nenv 
#%^!ifdef RANDOMEFFECTS
  random <- object$random
  if(!is.null(random))
    environment(random) <- nenv
#%^!endif  

  itags    <- object$Inter$itags
# The following are currently unused  
#  ninter   <- object$Inter$ninter
#  iused    <- object$Inter$iused
#  trivial  <- object$Inter$trivial

  # names of variables
  dat <- object$data
  datanames <- names(dat)
  dfnames <- summary(dat)$dfnames
  nondfnames <- datanames[!(datanames %in% dfnames)]
  
  # extract data-frame values
  dfdata <- as.data.frame(dat[, dfnames, drop=FALSE], warn=FALSE)
  
  # determine which interaction(s) are in force 
  answer <- impliedpresence(itags, iformula, dfdata, nondfnames)
#%^!ifdef RANDOMEFFECTS  
  if(!is.null(random)) {
    if("/" %in% all.names(random)) {
      ## hack since model.matrix doesn't handle "|" as desired
      rnd <- gsub("|", "/", pasteFormula(random), fixed=TRUE)
      random <- as.formula(rnd, env=environment(random))
    }
    answer2 <- impliedpresence(itags, random, dfdata, nondfnames)
    answer <- answer | answer2
  }
#%^!endif  
  colnames(answer) <- names(interaction)
  return(answer)
}

impliedcoefficients <- function(object, tag) {
  stopifnot(inherits(object, "mppm"))
  stopifnot(is.character(tag) && length(tag) == 1)
  fitobj      <- object$Fit$FIT
  Vnamelist   <- object$Fit$Vnamelist
# Not currently used:  
#  fitter      <- object$Fit$fitter
#  interaction <- object$Inter$interaction
#  ninteract   <- object$Inter$ninteract
#  trivial     <- object$Inter$trivial
#  iused       <- object$Inter$iused
  itags       <- object$Inter$itags
  if(!(tag %in% itags))
    stop(paste("Argument", dQuote("tag"),
               "is not one of the interaction names"))
  # (0) Set up
  # Identify the columns of the glm data frame
  # that are associated with this interpoint interaction 
  vnames <- Vnamelist[[tag]]
  if(!is.character(vnames))
    stop("Internal error - wrong format for vnames")
  # The answer is a matrix of coefficients,
  # with one row for each point pattern,
  # and one column for each vname
  answer <- matrix(, nrow=object$npat, ncol=length(vnames))
  colnames(answer) <- vnames

  # (1) make a data frame of covariates
  # Names of all columns in glm data frame
  allnames <- names(object$Fit$moadf)
  # Extract the design covariates
  df <- as.data.frame(object$data, warn=FALSE)
  # Names of all covariates other than design covariates
  othernames <- allnames[!(allnames %in% names(df))]
  # Add columns in which all other covariates are set to 0
  for(v in othernames) df[, v] <- 0
  
  # (2) evaluate linear predictor
  opt <- options(warn= -1)
  eta0 <- predict(fitobj, newdata=df, type="link")
  options(opt)
  
  # (3) for each vname in turn,
  # set the value of the vname to 1 and predict again
  for(j in seq(vnames)) {
    df[[vnames[j] ]] <- 1
    opt <- options(warn= -1)
    etaj <- predict(fitobj, newdata=df, type="link")
    options(opt)
    answer[ ,j] <- etaj - eta0
  }
  return(answer)
}



illegal.iformula <- local({

  illegal.iformula <- function(ifmla, itags, dfvarnames) {
    ## THIS IS TOO STRINGENT!
    ## Check validity of the interaction formula.
    ##  ifmla is the formula.
    ##  itags is the character vector of interaction names.
    ## Check whether the occurrences of `itags' in `iformula' are valid:
    ## e.g. no functions applied to `itags[i]'.
    ## Returns NULL if legal, otherwise a character string 
    stopifnot(inherits(ifmla, "formula"))
    stopifnot(is.character(itags))
    ## formula must not have a LHS
    if(length(ifmla) > 2)
      return("iformula must not have a left-hand side")
    ## variables in formula must be interaction tags or data frame variables
    varsinf <- variablesinformula(ifmla)
    if(!all(ok <- varsinf %in% c(itags, dfvarnames))) 
      return(paste(
                   ngettext(sum(!ok), "variable", "variables"),
                   paste(dQuote(varsinf[!ok]), collapse=", "),
                   "not allowed in iformula"))
    ## if formula uses no interaction tags, it's trivial
    if(!any(itags %in% variablesinformula(ifmla)))
      return(NULL)
    ## create terms object
    tt <- attributes(terms(ifmla))
    ## extract all variables appearing in the formula
    vars <- as.list(tt$variables)[-1]
    ##  nvars <- length(vars)
    varexprs <- lapply(vars, as.expression)
    varstrings <- sapply(varexprs, paste)
    ## Each variable may be a name or an expression
    v.is.name <- sapply(vars, is.name)
    ## a term may be an expression like sin(x), poly(x,y,degree=2)
    v.args <- lapply(varexprs, all.vars)
    ##  v.n.args <- sapply(v.args, length)
    v.has.itag <- sapply(lapply(v.args, "%in%", x=itags), any)
    ## interaction tags may only appear as names, not in functions
    if(any(nbg <- v.has.itag & !v.is.name))
      return(paste("interaction tags may not appear inside a function:",
                   paste(dQuote(varstrings[nbg]), collapse=", ")))
    ## Interaction between two itags is not defined
    ## Inspect the higher-order terms
    fax <- tt$factors
    if(prod(dim(fax)) == 0)
      return(NULL)
    ## rows are first order terms, columns are terms of order >= 1
    fvars <- rownames(fax)
    fterms <- colnames(fax)
    fv.args <- lapply(fvars, variablesintext)
    ft.args <- lapply(fterms, variables.in.term, 
                      factors=fax, varnamelist=fv.args)
    ft.itags <- lapply(ft.args, intersect, y=itags)
    if(any(sapply(ft.itags, length) > 1))
      return("Interaction between itags is not defined")
    return(NULL)
  }

  variables.in.term <- function(term, factors, varnamelist) {
    basis <- (factors[, term] != 0)
    unlist(varnamelist[basis])
  }
  
  illegal.iformula
})
