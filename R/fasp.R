#
#	fasp.R
#
#	$Revision: 1.30 $	$Date: 2013/01/24 07:47:34 $
#
#
#-----------------------------------------------------------------------------
#

# creator
fasp <- function(fns, which, formulae=NULL,
                 dataname=NULL, title=NULL, rowNames=NULL, colNames=NULL,
                 checkfv=TRUE) {
  stopifnot(is.list(fns))
  stopifnot(is.matrix(which))
  stopifnot(length(fns) == length(which))
  n   <- length(which)

  if(checkfv)
    for(i in seq_len(n))
      if(!is.fv(fns[[i]]))
        stop(paste("fns[[", i, "]] is not an fv object", sep=""))

  # set row and column labels
  if(!is.null(rowNames))
    rownames(which) <- rowNames
  if(!is.null(colNames))
    colnames(which) <- colNames

  if(!is.null(formulae)) {
    # verify format and convert to character vector
    formulae <- FormatFaspFormulae(formulae, "formulae")
    # ensure length matches length of "fns"
    if(length(formulae) == 1 && n > 1)
        # single formula - replicate it
        formulae <- rep(formulae, n)
    else 
        stopifnot(length(formulae) == length(which))
  }

  rslt <- list(fns=fns, 
               which=which, default.formula=formulae,
               dataname=dataname, title=title)
  class(rslt) <- "fasp"
  return(rslt)
}

# subset extraction operator

"[.fasp" <-
  function(x, I, J, drop=TRUE, ...) {

        verifyclass(x, "fasp")
        
        m <- nrow(x$which)
        n <- ncol(x$which)
        
        if(missing(I)) I <- 1:m
        if(missing(J)) J <- 1:n
        if(!is.vector(I) || !is.vector(J))
          stop("Subset operator is only implemented for vector indices")

        # determine index subset for lists 'fns', 'titles' etc
        included <- rep(FALSE, length(x$fns))
        w <- as.vector(x$which[I,J])
        if(length(w) == 0)
          stop("result is empty")
        included[w] <- TRUE

        # if only one cell selected, and drop=TRUE:
        if((sum(included) == 1) && drop)
          return(x$fns[included][[1]])
        
        # determine positions in shortened lists
        whichIJ <- x$which[I,J,drop=FALSE]
        newk <- cumsum(included)
        newwhich <- matrix(newk[whichIJ],
                           ncol=ncol(whichIJ), nrow=nrow(whichIJ))
        rownames(newwhich) <- rownames(x$which)[I]
        colnames(newwhich) <- colnames(x$which)[J]

        # default plotting formulae - could be NULL
        deform <- x$default.formula
        
        # create new fasp object
        Y <- fasp(fns      = x$fns[included],
                  formulae = if(!is.null(deform)) deform[included] else NULL,
                  which    = newwhich,
                  dataname = x$dataname,
                  title    = x$title)
        return(Y)
}

dim.fasp <- function(x) { dim(x$which) }

# print method

print.fasp <- function(x, ...) {
  verifyclass(x, "fasp")
  cat(paste("Function array (class", sQuote("fasp"), ")\n"))
  dim <- dim(x$which)
  cat(paste("Dimensions: ", dim[1], "x", dim[2], "\n"))
  cat(paste("Title:", if(is.null(x$title)) "(None)" else x$title, "\n"))
  invisible(NULL)
}

# other methods

dimnames.fasp <- function(x) {
  return(dimnames(x$which))
}

"dimnames<-.fasp" <- function(x, value) {
  w <- x$which
  dimnames(w) <- value
  x$which <- w
  return(x)
}

pool.fasp <- function(...) {
  Alist <- list(...)
  Yname <- short.deparse(sys.call())
  if(nchar(Yname) > 60) Yname <- paste(substr(Yname, 1, 40), "[..]")
  nA <-  length(Alist)
  if(nA == 0) return(NULL)
  # validate....
  # All arguments must be fasp objects
  notfasp <- !unlist(lapply(Alist, inherits, what="fasp"))
  if(any(notfasp)) {
    n <- sum(notfasp)
    why <- paste(ngettext(n, "Argument", "Arguments"),
                 commasep(which(notfasp)),
                 ngettext(n, "does not", "do not"),
                 "belong to the class",
                 dQuote("fasp"))
    stop(why)
  }
  # All arguments must have envelopes
  has.env <- function(z) {
    all(unlist(lapply(z$fns, inherits, what="envelope")))
  }
  notenv <- !unlist(lapply(Alist, has.env))
  if(any(notenv)) {
    n <- sum(notenv)
    why <- paste(ngettext(n, "Argument", "Arguments"),
                 commasep(which(notenv)),
                 ngettext(n, "does not", "do not"),
                 "contain envelope data")
    stop(why)
  }
  
  if(nA == 1) return(Alist[[1]])
  
  # All arguments must have the same dimensions
  witches <- lapply(Alist, function(z) { z$which })
  witch1 <- witches[[1]]
  same <- unlist(lapply(witches, identical, y=witch1))
  if(!all(same))
    stop("Function arrays do not have the same array dimensions")
  
  # OK.
  # Pool envelopes at each position
  result <- Alist[[1]]
  fns <- result$fns
  for(k in seq_along(fns)) {
    funks <- lapply(Alist, function(z, k) { z$fns[[k]] }, k=k)
    fnk <- do.call("pool.envelope", funks)
    attr(fnk, "einfo")$Yname <- Yname
    fns[[k]] <- fnk
  }
  result$fns <- fns

  return(result)
}

# other functions

FormatFaspFormulae <- function(f, argname) {
  # f should be a single formula object, a list of formula objects,
  # a character vector, or a list containing formulae and strings.
  # It will be converted to a character vector.
  
  zapit <- function(x, argname) {
    if(inherits(x, "formula")) deparse(x)
    else if(is.character(x)) x
    else stop(paste("The entries of",
                    sQuote(argname),
                    "must be formula objects or strings"))
  }

  result <-
    if(is.character(f))
      f
    else if(inherits(f, "formula"))
      deparse(f)
    else if(is.list(f))
      unlist(lapply(f, zapit, argname=argname))
    else stop(paste(sQuote(argname),
                    "should be a formula, a list of formulae,",
                    "or a character vector"))

  return(result)
}
