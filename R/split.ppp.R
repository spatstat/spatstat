#
# split.ppp.R
#
# $Revision: 1.25 $ $Date: 2015/02/01 01:24:58 $
#
# split.ppp and "split<-.ppp"
#
#########################################

split.ppp <- function(x, f = marks(x), drop=FALSE, un=NULL, ...) {
  verifyclass(x, "ppp")
  mf <- markformat(x)
  
  if(is.null(un))
    un <- missing(f) && (mf != "dataframe")

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
           dataframe={
             f <- fsplit <- firstfactor(marks(x))
             if(is.null(f))
               stop("Data frame of marks contains no factors")
           })
    splittype <- "factor"
  } else{
    # f was given
    fsplit <- f
    if(is.factor(f)) {
      splittype <- "factor"
    } else if(is.tess(f)) {
      # f is a tessellation: determine the grouping
      f <- marks(cut(x, fsplit))
      splittype <- "tess"
    } else if(is.owin(f)) {
      # f is a window: coerce to a tessellation
      W <- as.owin(x)
      fsplit <- tess(tiles=list(fsplit, setminus.owin(W, fsplit)),
                     window=W)
      f <- marks(cut(x, fsplit))
      splittype <- "tess"
    } else if(is.im(f)) {
      # f is an image: coerce to a tessellation
      fsplit <- tess(image=f)
      f <- marks(cut(x, fsplit))
      splittype <- "tess"
    } else if(is.character(f) && length(f) == 1) {
      # f is the name of a column of marks
      marx <- marks(x)
      if(is.data.frame(marx) && (f %in% names(marx))) 
        fsplit <- f <- marx[[f]]
      else
        stop(paste("The name", sQuote(f), "does not match any column of marks"))
      splittype <- "factor"
    } else 
      stop(paste("f must be",
                 "a factor, a tessellation, a window, an image,",
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
           tess = {
             # remove tiles that don't contain points
             fsplit <- fsplit[retain]
           },
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
  if(splittype == "tess") {
    til <- tiles(fsplit)
    for(i in seq_along(out))
      out[[i]]$window <- til[[i]]
  }
  class(out) <- c("splitppp", class(out))
  attr(out, "fsplit") <- fsplit
  attr(out, "fgroup") <- f
  return(out)
}

"split<-.ppp" <- function(x, f=marks(x), drop=FALSE, un=missing(f), 
                          ..., value) {
  verifyclass(x, "ppp")
  W <- x$window
  mf <- markformat(x)
  # evaluate `un' before assigning value of 'f'
  force(un)

  # validate assignment value
  stopifnot(is.list(value))
  if(!all(unlist(lapply(value, is.ppp))))
    stop(paste("Each entry of", sQuote("value"),
               "must be a point pattern"))

  ismark <- unlist(lapply(value, is.marked))
  if(any(ismark) && !all(ismark))
    stop(paste("Some entries of",
               sQuote("value"),
               "are marked, and others are unmarked"))
  vmarked <- all(ismark)

  # determine type of splitting
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
           dataframe={
             f <- fsplit <- firstfactor(marks(x))
             if(is.null(f))
               stop("Data frame of marks contains no factors")
           })
  } else {
    # f given
    fsplit <- f
    if(is.tess(f)) {
      # f is a tessellation: determine the grouping
      f <- marks(cut(x, fsplit))
    } else if(is.im(f)) {
      # f is an image: determine the grouping
      fsplit <- tess(image=f)
      f <- marks(cut(x, fsplit))
    } else if(is.character(f) && length(f) == 1) {
      # f is the name of a column of marks
      marx <- marks(x)
      if(is.data.frame(marx) && (f %in% names(marx))) 
        fsplit <- f <- marx[[f]]
      else
        stop(paste("The name", sQuote(f), "does not match any column of marks"))
    } else if(!is.factor(f))
      stop(paste("f must be",
                 "a factor, a tessellation, an image,",
                 "or the name of a column of marks"))
    if(length(f) != x$n)
      stop("length(f) must equal the number of points in x")
  } 
  #
  all.levels <- lev <- levels(f)
  if(!drop) 
    levtype <- "levels of f"
  else {
    levtype <- "levels which f actually takes"
    # remove components that don't contain points
    lev <- lev[table(f) > 0]
  }
  if(length(value) != length(lev))
      stop(paste("length of", sQuote("value"),
                 "should equal the number of",
                 levtype))

  # ensure value[[i]] is associated with lev[i]
  if(!is.null(names(value))) {
    if(!all(names(value) %in% as.character(lev)))
      stop(paste("names of", sQuote("value"), "should be levels of f"))
    value <- value[lev]
  }
  names(value) <- NULL

  # restore the marks, if they were discarded
  if(un && is.marked(x)) {
    if(vmarked)
      warning(paste(sQuote("value"), "contains marked point patterns:",
                    "this is inconsistent with un=TRUE; marks ignored."))
    for(i in seq_along(value)) 
      value[[i]] <- value[[i]] %mark% factor(lev[i], levels=all.levels)
  }

  # handle NA's in splitting factor
  if(any(isNA <- is.na(f))) {
    xNA <- x[isNA]
    if(un && is.marked(x)) 
      xNA <- xNA %mark% factor(NA, levels=all.levels)
    value <- append(value, list(xNA))
  }

  # put Humpty together again
  if(npoints(x) == length(f) &&
     length(levels(f)) == length(value) &&
     all(table(f) == sapply(value, npoints))) {
    ## exact correspondence
    out <- x
    for(i in seq_along(levels(f)))
      out[ f == lev[i] ] <- value[[i]]
  } else {
    out <- do.call(superimpose,c(value,list(W=W)))
  }
  return(out)
}


print.splitppp <- function(x, ...) {
  f <- attr(x, "fsplit")
  what <- if(is.tess(f)) "tessellation" else
          if(is.factor(f)) "factor" else "unknown data"
  cat(paste("Point pattern split by", what, "\n"))
  nam <- names(x)
  for(i in seq_along(x)) {
    cat(paste("\n", nam[i], ":\n", sep=""))
    print(x[[i]])
  }
  return(invisible(NULL))
}

summary.splitppp <- function(object, ...) {
  x <- lapply(object, summary, ...)
  class(x) <- "summary.splitppp"
  x
}

print.summary.splitppp <- function(x, ...) {
  class(x) <- "listof"
  print(x)
  invisible(NULL)
}

"[.splitppp" <- function(x, ...) {
  f <- attr(x, "fsplit")
  # invoke list method on x
  class(x) <- "list"
  y <- x[...]
  # then make it a 'splitppp' object too
  class(y) <- c("splitppp", class(y))
  if(is.tess(f)) {
    fsplit <- f[...]
  } else if(is.factor(f)) {
    lev <- levels(f)
    sublev <- lev[...]
    subf <- f[f %in% sublev]
    fsplit <- factor(subf, levels=lev)
  } else stop("Unknown splitting type")
  attr(y, "fsplit") <- fsplit
  y
}

"[<-.splitppp" <- function(x, ..., value) {
  if(!all(unlist(lapply(value, is.ppp))))
    stop("replacement value must be a list of point patterns")
  f <- attr(x, "fsplit")
  # invoke list method
  class(x) <- "list"
  x[...] <- value
  # then make it a 'splitppp' object too
  class(x) <- c("splitppp", class(x))
  if(is.tess(f)) {
    fsplit <- f
  } else if(is.factor(f)) {
    lev <- levels(f)
    fsplit <- factor(rep.int(lev, unlist(lapply(x, npoints))), levels=lev)
  }
  attr(x, "fsplit") <- fsplit
  x
}
  
density.splitppp <- function(x, ...) {
  as.listof(lapply(x, density, ...))
}

plot.splitppp <- function(x, ..., main) {
  if(missing(main)) main <- short.deparse(substitute(x))
  do.call("plot.listof",
          resolve.defaults(list(x=x, main=main),
                           list(...),
                           list(equal.scales=TRUE)))
}

as.layered.splitppp <- function(X) { do.call("layered", X) }

