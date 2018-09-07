#
# split.ppx.R
#
# $Revision: 1.6 $ $Date: 2018/09/07 05:49:15 $
#
# split.ppx etc
#
#########################################

split.ppx <- function(x, f = marks(x), drop=FALSE, un=NULL, ...) {
  stopifnot(inherits(x, "ppx"))
  mf <- markformat(x)
  if(is.null(un))
    un <- missing(f) && !(mf %in% c("dataframe", "hyperframe"))

  if(missing(f)) {
    # f defaults to marks of x
    switch(mf,
           none={
             stop("f is missing and there are no marks")
           },
           vector={
             if(!is.multitype(x)) 
               stop("f is missing and the pattern is not multitype")
             f <- fsplit <- marks(x)
           },
           hyperframe=,
           dataframe={
             f <- fsplit <- firstfactor(marks(x))
             if(is.null(f))
               stop("Marks do not include a factor")
           })
    splittype <- "factor"
  } else{
    # f was given
    fsplit <- f
    if(is.factor(f)) {
      splittype <- "factor"
    } else if(is.logical(f)) {
      splittype <- "factor"
      f <- factor(f)
    } else if(is.character(f) && length(f) == 1) {
      # f is the name of a column of marks
      marx <- marks(x)
      if((is.data.frame(marx) || is.hyperframe(marx))
         && (f %in% names(marx))) {
        fsplit <- f <- as.factor(marx[ ,f,drop=TRUE])
      } else
        stop(paste("The name", sQuote(f), "does not match any column of marks"))
      splittype <- "factor"
    } else 
      stop(paste("f must be",
                 "a factor,",
                 "or the name of a column of marks"))
    if(length(f) != npoints(x))
      stop("length(f) must equal the number of points in x")
  }

  # At this point
  # 'f' is a factor that can be used to separate the points
  # 'fsplit' is the object (either a factor or a tessellation)
  # that determines the split (and can be "un-split")

  lev <- levels(f)
  if(drop) {
    # remove components that don't contain points
    retain <- (table(f) > 0)
    lev <- lev[retain]
    switch(splittype,
           factor = {
             # delete levels that don't occur
             fsplit <- factor(fsplit, levels=lev)
           },
           stop("Internal error: wrong format for fsplit"))
  }

  # split the data
  out <- list()
  for(l in lev) 
    out[[paste(l)]] <- x[!is.na(f) & (f == l)]
  
  if(un)
     out <- lapply(out, unmark)
  class(out) <- c("splitppx", "anylist", class(out))
  attr(out, "fsplit") <- fsplit
  return(out)
}

print.splitppx <- function(x, ...) {
  f <- attr(x, "fsplit")
  what <- if(is.factor(f)) "factor" else "unknown data"
  cat(paste("Multidimensional point pattern split by", what, "\n"))
  nam <- names(x)
  for(i in seq_along(x)) {
    cat(paste("\n", nam[i], ":\n", sep=""))
    print(x[[i]])
  }
  return(invisible(NULL))
}

summary.splitppx <- function(object, ...) {
  x <- lapply(object, summary, ...)
  class(x) <- "summary.splitppx"
  x
}

print.summary.splitppx <- function(x, ...) {
  class(x) <- "anylist"
  print(x)
  invisible(NULL)
}

"[.splitppx" <- function(x, ...) {
  f <- attr(x, "fsplit")
  # invoke list method on x
  class(x) <- "list"
  y <- x[...]
  # then make it a 'splitppx' object too
  class(y) <- c("splitppx", class(y))
  if(is.factor(f)) {
    lev <- levels(f)
    sublev <- lev[...]
    subf <- f[f %in% sublev]
    fsplit <- factor(subf, levels=lev)
  } else stop("Unknown splitting type")
  attr(y, "fsplit") <- fsplit
  y
}

"[<-.splitppx" <- function(x, ..., value) {
  if(!all(unlist(lapply(value, is.ppx))))
    stop("replacement value must be a list of point patterns (ppx)")
  f <- attr(x, "fsplit")
  # invoke list method
  class(x) <- "list"
  x[...] <- value
  # then make it a 'splitppx' object too
  class(x) <- c("splitppx", class(x))
  if(is.factor(f)) {
    lev <- levels(f)
    fsplit <- factor(rep.int(lev, unlist(lapply(x, npoints))), levels=lev)
  }
  attr(x, "fsplit") <- fsplit
  x
}
  
