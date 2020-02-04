#
# Determine which 'canonical variables' depend on a supplied covariate
#
#   $Revision: 1.9 $  $Date: 2020/02/04 03:26:37 $
#

model.depends <- function(object) {
  # supplied covariates
  fo <- formula(object)
  if(length(as.list(fo)) == 3) {
    # formula has a response: strip it
    fo <- fo[-2]
  }
  covars <- variablesinformula(fo)
  # canonical covariates 
  mm <- model.matrix(object)
  ass <- attr(mm, "assign")
  # model terms
  tt <- terms(object)
  lab <- attr(tt, "term.labels")
  # 'ass' maps canonical covariates to 'lab'
  # determine which canonical covariate depends on which supplied covariate
  depends <- matrix(FALSE, length(ass), length(covars))
  for(i in seq(along=ass)) {
    if(ass[i] == 0) # 0 is the intercept term
      depends[i,] <- FALSE
    else {
      turm <- lab[ass[i]]
      depends[i, ] <- covars %in% all.vars(parse(text=turm))
    }
  }
  rownames(depends) <- colnames(mm)
  colnames(depends) <- covars
  # detect offsets
  if(!is.null(oo <- attr(tt, "offset")) && ((noo <- length(oo)) > 0)) {
    # entries of 'oo' index the list of variables in terms object
    vv <- attr(tt, "variables")
    offdep <- matrix(FALSE, noo, length(covars))
    offnms <- character(noo)
    for(i in seq_len(noo)) {
      offseti <- languageEl(vv, oo[i] + 1)
      offdep[i, ] <- covars %in% all.vars(offseti)
      offnms[i] <- deparse(offseti)
    }
    rownames(offdep) <- offnms
    colnames(offdep) <- covars
    attr(depends, "offset") <- offdep
  }
  return(depends)
}

model.is.additive <- function(object) {
  dep <- model.depends(object)
  hit <- t(dep) %*% dep
  diag(hit) <- 0
  ok <- all(hit == 0)
  return(ok)
}

model.covariates <- function(object, fitted=TRUE, offset=TRUE) {
  md <- model.depends(object)
  nm <- colnames(md)
  keep <- rep.int(FALSE, length(nm))
  # variables used in formula with coefficients
  if(fitted) keep <- apply(md, 2, any)
  # variables used in offset
  if(offset) {
    oo <- attr(md, "offset")
    if(!is.null(oo)) 
      keep <- keep | apply(oo, 2, any)
  }
  return(nm[keep])
}

has.offset.term <- function(object) {
  # model terms
  tt <- terms(object)
  oo  <- attr(tt, "offset")
  return(!is.null(oo) && (length(oo) > 0))
}

has.offset <- function(object) {
  has.offset.term(object) || !is.null(model.offset(model.frame(object)))
}

check.separable <- function(dmat, covname, isconstant, fatal=TRUE) {
  #' Determine whether the effect of 'covname' is separable from other terms.
  #' dmat = model.depends(model)
  #' Find covariates entangled with 'covname' in the model
  relevant <- dmat[, covname]
  othercov <- (colnames(dmat) != covname)
  conflict <- dmat[relevant, othercov, drop=FALSE]
  if(!any(conflict)) return(TRUE)
  #' names of entangled covariates
  entangled <- colnames(conflict)[matcolany(conflict)]
  #' not problematic if constant
  if(is.null(names(isconstant))) names(isconstant) <- colnames(dmat)
  ok <- unlist(isconstant[entangled])
  conflict[ , ok] <- FALSE
  if(!any(conflict)) return(TRUE)
  #' there are conflicts
  conflictterms <- matrowany(conflict)
  conflictcovs  <- matcolany(conflict)
  whinge <- paste("The covariate", sQuote(covname),
                  "cannot be separated from the",
                  ngettext(sum(conflictcovs), "covariate", "covariates"),
                  commasep(sQuote(colnames(conflict)[conflictcovs])),
                  "in the model",
                  ngettext(sum(conflictterms), "term", "terms"),
                  commasep(sQuote(rownames(conflict)[conflictterms])))
  if(fatal) stop(whinge, call.=FALSE)
  warning(whinge, call.=FALSE)
  return(FALSE)
}
