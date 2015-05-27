## hierarchy.R
##
## Support functions for hierarchical interactions
##
##  $Revision: 1.1 $  $Date: 2015/05/26 08:39:56 $

hierarchicalordering <- function(i, s) {
  s <- as.character(s)
  if(inherits(i, "hierarchicalordering")) {
    ## already a hierarchical ordering
    if(length(s) != length(i$labels))
      stop("Tried to change the number of types in the hierarchical order")
    i$labels <- s
    return(i)
  }
  n <- length(s)
  possible <- if(is.character(i)) s else seq_len(n)
  j <- match(i, possible)
  if(any(uhoh <- is.na(j)))
    stop(paste("Unrecognised",
               ngettext(sum(uhoh), "level", "levels"),
               commasep(sQuote(i[uhoh])),
               "amongst possible levels",
               commasep(sQuote(s))))
  if(length(j) < n)
    stop("Ordering is incomplete")
  ord <- order(j)
  m <- matrix(, n, n)
  rel <- matrix(ord[row(m)] <= ord[col(m)], n, n)
  dimnames(rel) <- list(s, s)
  x <- list(indices=j, ordering=ord, labels=s, relation=rel)
  class(x) <- "hierarchicalordering"
  x
}

print.hierarchicalordering <- function(x, ...) {
  splat(x$labels[x$indices], collapse=" ~> ")
  invisible(NULL)
}
                     
hiermat <- function (x, h) 
{
  stopifnot(is.matrix(x))
  isna <- is.na(x)
  x[] <- as.character(x)
  x[isna] <- ""
  if(inherits(h, "hierarchicalordering")) ## allows h to be NULL, etc
    x[!(h$relation)] <- ""
  return(noquote(x))
}
