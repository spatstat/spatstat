#
#  unnormdensity.R
#
#  $Revision: 1.3 $  $Date: 2011/10/22 04:47:10 $
#

unnormdensity <- function(x, ..., weights=NULL) {
  if(any(!nzchar(names(list(...)))))
    stop("All arguments must be named (tag=value)")
  if(is.null(weights)) {
    out <- do.call.matched("density.default", c(list(x=x), list(...)))
  } else if(all(weights == 0)) {
    # result is zero
    out <- do.call.matched("density.default", c(list(x=x), list(...)))
    out$y <- 0 * out$y
  } else if(all(weights >= 0)) {
    # all masses are nonnegative
    w <- weights
    totmass <- sum(w)
    out <- do.call.matched("density.default",
                           c(list(x=x),
                             list(...),
                             list(weights=w/totmass)))
    out$y <- out$y * totmass
  } else if(all(weights <= 0)) {
    # all masses are nonpositive
    w <- (- weights)
    totmass <- sum(w)
    out <- do.call.matched("density.default",
                           c(list(x=x),
                             list(...),
                             list(weights=w/totmass)))
    out$y <- out$y * (- totmass)
  } else {
    # mixture of positive and negative masses
    w <- weights
    wabs <- abs(w)
    wpos <- pmax(0, w)
    wneg <- - pmin(0, w)
    # determine bandwidth using absolute masses
    dabs <- do.call.matched("density.default",
                            c(list(x=x),
                              list(...),
                              list(weights=wabs/sum(wabs))))
    bw <- dabs$bw
    # compute densities for positive and negative masses separately
    outpos <- do.call.matched("density.default",
                      resolve.defaults(list(x=x),
                                       list(bw=bw, adjust=1),
                                       list(weights=wpos/sum(wpos)),
                                       list(...),
                                       .StripNull=TRUE))
    outneg <- do.call.matched("density.default",
                      resolve.defaults(list(x=x),
                                       list(bw=bw, adjust=1),
                                       list(weights=wneg/sum(wneg)),
                                       list(...),
                                       .StripNull=TRUE))
    # combine
    out <- outpos
    out$y <- sum(wpos) * outpos$y - sum(wneg) * outneg$y
  }
  out$call <- match.call()
  return(out)
}

