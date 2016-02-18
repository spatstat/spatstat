#'
#'   factors.R
#'
#'  Tools for manipulating factors and factor-valued things
#'
#'  $Revision: 1.3 $  $Date: 2015/03/28 10:22:42 $

relevel.im <- function(x, ref, ...) {
  if(x$type != "factor")
    stop("Only valid for factor-valued images")
  x[] <- relevel(x[], ref, ...)
  return(x)
}

relevel.ppp <- relevel.ppx <- function(x, ref, ...) {
  stopifnot(is.multitype(x))
  marks(x) <- relevel(marks(x), ref, ...)
  return(x)
}

mergeLevels <- function(.f, ...) {
  if(is.im(.f)) {
    aa <- mergeLevels(.f[], ...)
    .f[] <- aa
    return(.f)
  }
  if(is.multitype(.f)) {
    marks(.f) <- mergeLevels(marks(.f), ...)
    return(.f)
  }
  stopifnot(is.factor(.f))
  map <- list(...)
  n <- length(map)
  if(n == 0) return(.f)
  # mapping for 'other'
  if(any(isnul <- (lengths(map) == 0))) {
    if(sum(isnul) > 1)
      stop("At most one argument should be NULL or character(0)")
    otherlevels <- setdiff(levels(.f), unlist(map))
    map[[which(isnul)]] <- otherlevels
  }
  newlevels <- names(map)
  oldlevels <- levels(.f)
  mappedlevels <- unlist(map)
  if(sum(nzchar(newlevels)) != n)
    stop("Arguments must be in the form name=value")
  if(!all(mappedlevels %in% oldlevels))
    stop("Argument values must be levels of .f")
  ## construct mapping
  fullmap <- oldlevels
  for(i in seq_len(n)) {
    relevant <- oldlevels %in% map[[i]]
    fullmap[relevant] <- newlevels[i]
  }
  ## apply mapping
  newf <- factor(fullmap[.f], levels=unique(fullmap))
  return(newf)
}



