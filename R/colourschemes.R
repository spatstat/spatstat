#
#  colourschemes.R
#
#  $Revision: 1.4 $  $Date: 2019/01/29 04:32:44 $
#

beachcolourmap <- function(range, ...) {
  col <- beachcolours(range, ...)
  z <- colourmap(col, range=range)
  return(z)
}

beachcolours <- function(range, sealevel = 0, monochrome=FALSE,
                         ncolours=if(monochrome) 16 else 64,
                         nbeach=1) {
  check.range(range)
  stopifnot(all(is.finite(range)))
  check.1.real(sealevel)
  check.in.range(sealevel, range)
  check.1.integer(ncolours)
  stopifnot(ncolours >= 3)
  check.1.integer(nbeach)
  stopifnot(nbeach >= 0)
  stopifnot(nbeach <= ncolours + 2)
  if(monochrome)
    return(grey(seq(from=0,to=1,length.out=ncolours)))
  depths <- range[1L]
  peaks <- range[2L]
  dv <- diff(range)/(ncolours - 1L)
  epsilon <- nbeach * dv/2
  lowtide <- max(sealevel - epsilon, depths)
  hightide <-  min(sealevel + epsilon, peaks)
  nsea <- max(0L, floor((lowtide - depths)/dv))
  nland <- max(0L, floor((peaks  - hightide)/dv))
  discrep <- nsea + nland + nbeach - ncolours
  if(discrep != 0) {
    dd <- abs(discrep)
    ss <- as.integer(-sign(discrep))
    smallhalf <- dd/2L
    largehalf <- dd - smallhalf
    if(nsea < nland) {
      nsea <- nsea + ss * smallhalf
      nland <- nland + ss * largehalf
    } else {
      nland <- nland + ss * smallhalf
      nsea <- nsea + ss * largehalf
    }      
    if(nsea + nland + nbeach != ncolours)
      warning("Internal error: incorrect adjustment of length in beachcolours")
  }
  colours <- character(0)
  if(nsea > 0)  colours <- rev(rainbow(nsea, start=3/6,end=4/6)) # cyan/blue
  if(nbeach > 0)  colours <- c(colours,
                               rev(rainbow(nbeach, start=3/12,end=5/12))) # green
  if(nland > 0)  colours <- c(colours,
                              rev(rainbow(nland, start=0, end=1/6)))  # red/yellow
  return(colours)
}
