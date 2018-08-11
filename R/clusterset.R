#
#   clusterset.R
#
#   Allard-Fraley estimator of cluster region
#
#   $Revision: 1.12 $  $Date: 2016/02/16 01:39:12 $
#

clusterset <- function(X, what=c("marks", "domain"),
                       ...,
                       verbose=TRUE,
                       fast=FALSE,
                       exact=!fast) {
  stopifnot(is.ppp(X))
  what <- match.arg(what, several.ok=TRUE)
  if(!missing(exact)) stopifnot(is.logical(exact))
  if(fast && exact)
    stop("fast=TRUE is incompatible with exact=TRUE")
  # compute duplication exactly as in deldir, or the universe will explode
  X <- unique(unmark(X), rule="deldir", warn=TRUE)
  n <- npoints(X)
  W <- as.owin(X)
  # discretised Dirichlet tessellation
  if(verbose) cat("Computing Dirichlet tessellation...")
  if(fast || !exact)
    cellid <- as.im(nnfun(X), ...)
  # compute tile areas
  if(fast) {
    a <- table(factor(as.vector(as.matrix(cellid)), levels=1:n))
    if(verbose) cat("done.\n")
    a <- a + 0.5
    A <- sum(a)
  } else {
    d <- dirichlet(X)
    if(verbose) cat("done.\n")
    D <- tiles(d)
    suppressWarnings(id <- as.integer(names(D)))
    if(anyNA(id) && ("marks" %in% what))
      stop("Unable to map Dirichlet tiles to data points")
    A <- area(W)
    a <- unlist(lapply(D, area))
  }
  # determine optimal selection of tiles
  ntile <- length(a)
  o <- order(a)
  b <- cumsum(a[o])
  m <- seq_len(ntile)
  logl <- -n * log(n) + m * log(m/b) + (n-m) * log((n-m)/(A-b))
  mopt <- which.max(logl)
  picked <- o[seq_len(mopt)]
  ## map tiles to points
  if(!fast) picked <- id[picked]
  ## logical vector
  is.picked <- rep.int(FALSE, n)
  is.picked[picked] <- TRUE
  # construct result
  out <- list(marks=NULL, domain=NULL)
  if("marks" %in% what) {
    ## label points
    yesno <- factor(ifelse(is.picked, "yes", "no"), levels=c("no", "yes"))
    out$marks <- X %mark% yesno
  }
  if("domain" %in% what) {
    if(verbose) cat("Computing cluster set...")
    if(exact) {
      domain <- do.call(union.owin, unname(D[is.picked]))
      domain <- rebound.owin(domain, as.rectangle(W))
    } else {
      domain <- eval.im(is.picked[cellid])
    }
    out$domain <- domain
    if(verbose) cat("done.\n")
  }
  out <- if(length(what) == 1L) out[[what]] else out
  return(out)
}
