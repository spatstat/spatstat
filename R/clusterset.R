#
#   clusterset.R
#
#   Allard-Fraley estimator of cluster region
#
#   $Revision: 1.4 $  $Date: 2012/10/11 10:30:52 $
#

clusterset <- function(X, result=c("marks", "domain"),
                       ...,
                       verbose=TRUE,
                       fast=FALSE,
                       exact=!fast && spatstat.options("gpclib")) {
  stopifnot(is.ppp(X))
  result <- match.arg(result)
  if(!missing(exact)) {
    stopifnot(is.logical(exact))
    if(exact && result == "domain" && !spatstat.options("gpclib"))
      stop(paste("Cannot perform exact calculation because gpclib is disabled;",
                 "see help(licence.polygons)"))
  }
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
    if(any(is.na(id)) && result == "marks")
      stop("Unable to map Dirichlet tiles to data points")
    A <- area.owin(W)
    a <- unlist(lapply(D, area.owin))
  }
  # determine optimal selection of tiles
  ntile <- length(a)
  o <- order(a)
  b <- cumsum(a[o])
  m <- seq_len(ntile)
  logl <- -n * log(n) + m * log(m/b) + (n-m) * log((n-m)/(A-b))
  mopt <- which.max(logl)
  picked <- o[seq_len(mopt)]
  # construct result
  switch(result,
         marks = {
           # map tiles to points
           if(!fast) picked <- id[picked]
           # label points
           is.picked <- rep("no", n)
           is.picked[picked] <- "yes"
           is.picked <- factor(is.picked, levels=c("no", "yes"))
           out <- X %mark% is.picked
         },
         domain = {
           if(exact) {
             out <- do.call("union.owin", unname(D[picked]))
           } else {
             is.picked <- rep(FALSE, n)
             is.picked[picked] <- TRUE
             out <- eval.im(is.picked[cellid])
           }
         })
  return(out)
}
