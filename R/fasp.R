#
#	fasp.R
#
#	$Revision: 1.36 $	$Date: 2020/11/24 01:38:24 $
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
    if(length(formulae) == 1L && n > 1L)
        # single formula - replicate it
        formulae <- rep.int(formulae, n)
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
        included <- rep.int(FALSE, length(x$fns))
        w <- as.vector(x$which[I,J])
        if(length(w) == 0)
          stop("result is empty")
        included[w] <- TRUE

        # if only one cell selected, and drop=TRUE:
        if((sum(included) == 1L) && drop)
          return(x$fns[included][[1L]])
        
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
  cat(paste("Dimensions: ", dim[1L], "x", dim[2L], "\n"))
  cat(paste("Title:", if(is.null(x$title)) "(None)" else x$title, "\n"))
  invisible(NULL)
}

# other methods

as.fv.fasp <- function(x) do.call(cbind.fv, x$fns)

dimnames.fasp <- function(x) {
  return(dimnames(x$which))
}

"dimnames<-.fasp" <- function(x, value) {
  w <- x$which
  dimnames(w) <- value
  x$which <- w
  return(x)
}


## other functions

FormatFaspFormulae <- local({

  zapit <- function(x, argname) {
    if(inherits(x, "formula")) deparse(x)
    else if(is.character(x)) x
    else stop(paste("The entries of",
                    sQuote(argname),
                    "must be formula objects or strings"))
  }

  FormatFaspFormulae <- function(f, argname) {
    ## f should be a single formula object, a list of formula objects,
    ## a character vector, or a list containing formulae and strings.
    ## It will be converted to a character vector.
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

  FormatFaspFormulae
})


