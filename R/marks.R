#
# marks.R
#
#   $Revision: 1.43 $   $Date: 2016/02/16 01:39:12 $
#
# stuff for handling marks
#
#

marks <- function(x, ...) {
  UseMethod("marks")
}

marks.default <- function(x, ...) { NULL }

# The 'dfok' switch is temporary
# while we convert the code to accept data frames of marks

marks.ppp <- function(x, ..., dfok=TRUE, drop=TRUE) {
  ma <- x$marks
  if((is.data.frame(ma) || is.matrix(ma))) {
    if(!dfok)
      stop("Sorry, not implemented when the marks are a data frame")
    if(drop && ncol(ma) == 1)
      ma <- ma[,1,drop=TRUE]
  }
  return(ma)
}

# ------------------------------------------------------------------

"marks<-" <- function(x, ..., value) {
  UseMethod("marks<-")
}

"marks<-.ppp" <- function(x, ..., dfok=TRUE, drop=TRUE, value) {
  np <- npoints(x)
  m <- value
  switch(markformat(m),
         none = {
           return(unmark(x))
         },
         vector = {
           # vector of marks
           if(length(m) == 1) m <- rep.int(m, np)
           else if(np == 0) m <- rep.int(m, 0) # ensures marked pattern obtained
           else if(length(m) != np) stop("number of points != number of marks")
           marx <- m
         },
         dataframe = {
           if(!dfok)
             stop("Sorry, data frames of marks are not yet implemented")
           m <- as.data.frame(m)
           # data frame of marks
           if(ncol(m) == 0) {
             # no mark variables
             marx <- NULL
           } else {
             # marks to be attached
             if(nrow(m) == np) {
               marx <- m
             } else {
               # lengths do not match
               if(nrow(m) == 1 || np == 0) {
               # replicate data frame
                 marx <- as.data.frame(lapply(as.list(m),
                                              function(x, k) { rep.int(x, k) },
                                              k=np))
               } else
               stop("number of rows of data frame != number of points")
             }
             # convert single-column data frame to vector?
             if(drop && ncol(marx) == 1)
               marx <- marx[,1,drop=TRUE]
           }
         },
         hyperframe = 
         stop("Hyperframes of marks are not supported in ppp objects; use ppx"),
         stop("Format of marks is not understood")
         )
  # attach/overwrite marks
  Y <- ppp(x$x,x$y,window=x$window,marks=marx, check=FALSE, drop=drop)
  return(Y)
}

"%mark%" <- setmarks <- function(x,value) {
  marks(x) <- value
  return(x)
}

# -------------------------------------------------

markformat <- function(x) {
  UseMethod("markformat")
}

markformat.ppp <- function(x) {
  mf <- x$markformat
  if(is.null(mf)) 
    mf <- markformat(marks(x))
  return(mf)
}

markformat.default <- function(x) {
  if(is.null(x)) return("none")
  if(is.null(dim(x))) {
    if(is.vector(x) || is.factor(x) || is.atomic(x)) return("vector")
    if(inherits(x, "POSIXt") || inherits(x, "Date")) return("vector")
  }
  if(is.data.frame(x) || is.matrix(x)) return("dataframe")
  if(is.hyperframe(x)) return("hyperframe")
  if(inherits(x, c("solist", "anylist", "listof"))) return("list")
  stop("Mark format not understood")
}

# ------------------------------------------------------------------

"is.marked" <-
function(X, ...) {
  UseMethod("is.marked")
}

"is.marked.ppp" <-
function(X, na.action="warn", ...) {
  marx <- marks(X, ...)
  if(is.null(marx))
    return(FALSE)
  if((length(marx) > 0) && anyNA(marx)) {
    gripe <- paste("some mark values are NA in the point pattern",
                   short.deparse(substitute(X)))
    switch(na.action,
           warn = warning(gripe, call.=FALSE),
           fatal = stop(gripe, call.=FALSE),
           ignore = {}
           )
  }
  return(TRUE)
}

"is.marked.default" <-
  function(...) { return(!is.null(marks(...))) }


# ------------------------------------------------------------------

is.multitype <- function(X, ...) {
  UseMethod("is.multitype")
}

is.multitype.default <- function(X, ...) {
  m <- marks(X)
  if(is.null(m))
    return(FALSE)
  if(!is.null(dim(m))) {
    # should have a single column
    if(dim(m)[2] != 1)
      return(FALSE)
    m <- m[,1,drop=TRUE]
  }
  return(is.factor(m))
}

is.multitype.ppp <- function(X, na.action="warn", ...) {
  marx <- marks(X, dfok=TRUE)
  if(is.null(marx))
    return(FALSE)
  if((is.data.frame(marx) || is.hyperframe(marx)) && ncol(marx) > 1)
    return(FALSE)
  if(!is.factor(marx))
    return(FALSE)
  if((length(marx) > 0) && anyNA(marx))
    switch(na.action,
           warn = {
             warning(paste("some mark values are NA in the point pattern",
                           short.deparse(substitute(X))))
           },
           fatal = {
             return(FALSE)
           },
           ignore = {}
           )
  return(TRUE)
}

# ------------------------------------------------------------------

unmark <- function(X) {
  UseMethod("unmark")
}

unmark.ppp <- function(X) {
  X$marks <- NULL
  X$markformat <- "none"
  return(X)
}

unmark.splitppp <- function(X) {
  Y <- lapply(X, unmark.ppp)
  class(Y) <- c("splitppp", class(Y))
  return(Y)
}

##### utility functions for subsetting & combining marks #########


marksubset <- function(x, index, format=NULL) {
  if(is.null(format)) format <- markformat(x)
  switch(format,
         none={return(NULL)},
         list=,
         vector={return(x[index])},
         hyperframe=,
         dataframe={return(x[index,,drop=FALSE])},
         stop("Internal error: unrecognised format of marks"))
}

"%msub%" <- marksubsetop <- function(x,i) { marksubset(x, i) }

"%mrep%" <- markreplicateop <- function(x,n) { 
  format <- markformat(x)
  switch(format,
         none={return(NULL)},
         list=,
         vector={ return(rep.int(x,n))},
         dataframe={
           return(as.data.frame(lapply(x, rep, times=n)))
         },
         hyperframe={
           xcols <- as.list(x)
           repxcols <- lapply(xcols, rep, times=n)
           return(do.call(hyperframe, repxcols))
         },
         stop("Internal error: unrecognised format of marks"))
}

"%mapp%" <- markappendop <- function(x,y) { 
  fx <- markformat(x)
  fy <- markformat(y)
  agree <- (fx == fy)
  if(all(c(fx,fy) %in% c("dataframe", "hyperframe")))
    agree <- agree && identical(names(x),names(y)) 
  if(!agree)
    stop("Attempted to concatenate marks that are not compatible")
  switch(fx,
         none   = { return(NULL) },
         vector = {
           if(is.factor(x) || is.factor(y))
             return(cat.factor(x,y))
           else return(c(x,y))
         },
         hyperframe=,
         dataframe = { return(rbind(x,y)) },
         list = {
           z <- append(x,y)
           z <- as.solist(z, demote=TRUE)
           return(z)
         },
         stop("Internal error: unrecognised format of marks"))
}

markappend <- function(...) {
  # combine marks from any number of patterns
  marxlist <- list(...)
  # check on compatibility of marks
  mkfmt <- sapply(marxlist,markformat)
  if(length(ufm <- unique(mkfmt))>1)
    stop(paste("Cannot append marks of different formats:",
               commasep(sQuote(ufm))),
         call.=FALSE)
  mkfmt <- mkfmt[1]
  # combine the marks
  switch(mkfmt,
         none = {
           return(NULL)
         },
         vector = {
           marxlist <- lapply(marxlist,
                              function(x){as.data.frame.vector(x,nm="v1")})
           marx <- do.call(rbind, marxlist)[,1]
           return(marx)
         },
         hyperframe =,
         dataframe = {
           # check compatibility of data frames
           # (this is redundant but gives more helpful message)
           nama <- lapply(marxlist, names)
           dims <- lengths(nama)
           if(length(unique(dims)) != 1)
             stop("Data frames of marks have different column dimensions.")
           samenames <- unlist(lapply(nama,
                                      function(x,y) { identical(x,y) },
                                      y=nama[[1]]))
           if(!all(samenames))
             stop("Data frames of marks have different names.\n")
           marx <- do.call(rbind, marxlist)
           return(marx)
         },
         list = {
           marx <- do.call(c, marxlist)
           marx <- as.solist(marx, demote=TRUE) 
           return(marx)
         })
  stop("Unrecognised mark format")
}

markcbind <- function(...) {
  # cbind several columns of marks
  marxlist <- list(...)
  mkfmt <- unlist(lapply(marxlist, markformat))
  if(any(vacuous <- (mkfmt == "none"))) {
    marxlist <- marxlist[!vacuous]
    mkfmt    <- mkfmt[!vacuous]
  }
  if(any(isvec <- (mkfmt == "vector"))) {
    ## convert vectors to data frames with invented names
    for(i in which(isvec)) {
      mi <- as.data.frame(marxlist[i])
      colnames(mi) <- paste0("V", i)
      marxlist[[i]] <- mi
    }
    mkfmt[isvec] <- "dataframe"
  }
  if(all(mkfmt == "dataframe")) {
    ## result is a data frame
    marx <- do.call(data.frame, marxlist)
  } else {
    ## result is a hyperframe
    if(!all(ishyp <- (mkfmt == "hyperframe"))) 
      marxlist[!ishyp] <- lapply(marxlist[!ishyp], as.hyperframe)
    marx <- do.call(hyperframe, marxlist)
  }
  return(marx)
}

# extract only the columns of (passably) numeric data from a data frame
numeric.columns <- function(M, logical=TRUE, others=c("discard", "na")) {
  others <- match.arg(others)
  M <- as.data.frame(M)
  if(ncol(M) == 1)
    colnames(M) <- NULL
  process <- function(z, logi, other) {
    if(is.numeric(z)) return(z)
    if(logi && is.logical(z)) return(as.integer(z))
    switch(other,
           na=rep.int(NA_real_, length(z)),
           discard=NULL,
           NULL)
  }
  Mprocessed <- lapply(M, process, logi=logical, other=others)
  isnul <- unlist(lapply(Mprocessed, is.null))
  if(all(isnul)) {
    # all columns have been removed
    # return a data frame with no columns
    return(as.data.frame(matrix(, nrow=nrow(M), ncol=0)))
  }
  Mout <- do.call(data.frame, Mprocessed[!isnul])
  if(ncol(M) == 1 && ncol(Mout) == 1)
    colnames(Mout) <- NULL
  return(Mout)
}

coerce.marks.numeric <- function(X, warn=TRUE) {
  marx <- marks(X)
  if(is.null(dim(marx))) {
    if(is.factor(marx)) {
      if(warn) warning("Factor-valued marks were converted to integer codes",
                       call.=FALSE)
      marx <- as.integer(marx)
      return(X %mark% marx)
    }
  } else {
    marx <- as.data.frame(marx)
    if(any(fax <- unlist(lapply(marx, is.factor)))) {
      if(warn) {
        nf <- sum(fax)
        whinge <- paste("Factor-valued mark",
                        ngettext(nf, "variable", "variables"),
                        commasep(sQuote(colnames(marx)[fax])),
                        ngettext(nf, "was", "were"),
                        "converted to integer codes")
        warning(whinge, call.=FALSE)
      }
      marx[fax] <- as.data.frame(lapply(marx[fax], as.integer))
      return(X %mark% marx)
    }
  }
  return(X)
}
