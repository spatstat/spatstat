#
#  reduceformula.R
#
#  $Revision: 1.3 $   $Date: 2007/04/02 06:28:17 $
#
# delete variable from formula 
#
#......................................................
#

reduceformula <- function(fmla, deletevar, verbose=FALSE) {
  # removes the variable `deletevar' from the formula `fmla'
  # returns a simplified formula, or NULL if it can't simplify.
  stopifnot(inherits(fmla, "formula"))
  stopifnot(is.character(deletevar) && length(deletevar) == 1)
  if(!(deletevar %in% all.vars(as.expression(fmla)))) {
    if(verbose)
      message(paste("The formula does not involve", dQuote(deletevar),
                    "and is therefore unchanged"))
    return(fmla)
  }
  lhs <- if(length(fmla) < 3) NULL else fmla[[2]]
  # create terms object
  tt <- attributes(terms(fmla))
  formula.has.intercept <- (tt$intercept == 1)
  # extract all variables appearing in the model
  vars <- as.list(tt$variables)[-1]
  nvars <- length(vars)
  varstrings <- sapply(vars, function(x) paste(as.expression(x)))
  # identify any offsets
  offs <- tt$offset
  v.is.offset <- if(!is.null(offs)) (1:nvars) %in% offs else rep(FALSE, nvars)
  # remove the response
  repo <- tt$response
  if(repo != 0) {
    vars <- vars[-repo]
    varstrings <- varstrings[-repo]
    v.is.offset <- v.is.offset[-repo]
  }
  # a term may be a variable name
  v.is.name <- sapply(vars, is.name)
  # a term may be an expression like sin(x), poly(x,y,degree=2)
  v.args <- lapply(vars, function(x) all.vars(as.expression(x)))
  v.has.delete <- sapply(v.args,
                         function(x,d) { d %in% x },
                         d=deletevar)
  v.has.other <- sapply(v.args,
                        function(x,d) { !all(x == d) },
                        d=deletevar)
  v.is.mixed <- v.has.delete & v.has.other
  # we can't handle mixed terms like sin(x-d), poly(x,d)
  # where d is to be deleted. Handling these would require
  # knowledge about the functions sin and poly.
  if(any(v.is.mixed)) {
    nmixed <- sum(v.is.mixed)
    if(verbose)
      message(paste("Don't know how to reduce the",
              ngettext(nmixed, "term", "terms"),
              paste(dQuote(varstrings[v.is.mixed]), collapse=",")))
    return(NULL)
  }
  # OK. We have identified all first order terms to be deleted.
  condemned.names <- varstrings[v.has.delete]
  # Determine the terms of all orders that include these first order terms
  # (1) terms with model coefficients
  fax <- tt$factors
  if(prod(dim(fax)) == 0)
    retained.terms <- character(0)
  else {
    # Rows are first order terms 
    condemned.row <- rownames(fax) %in% condemned.names
    # Columns are the terms of all orders
    allterms <- colnames(fax)
    # Find all columns containing a 1 in a row that is to be deleted
    if(any(condemned.row)) {
      condemned.column <- apply(fax[condemned.row, , drop=FALSE] != 0, 2, any)
      retained.terms <- allterms[!condemned.column]
    } else retained.terms <- allterms
  }
  # (2) offsets if any
  if(any(v.is.offset))
    retained.terms <- c(retained.terms,
                        varstrings[v.is.offset & !v.has.delete])
  # (3) intercept forced?
  if(length(retained.terms) == 0)
    retained.terms <- "1"
  
  # OK. Cut-and-paste
  f <- paste(lhs, "~", paste(retained.terms, collapse=" + "))
  return(as.formula(f))
} 

