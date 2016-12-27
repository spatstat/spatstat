#
#  colourschemes.R
#
#  $Revision: 1.3 $  $Date: 2013/07/17 04:53:48 $
#

beachcolourmap <- function(range, ...) {
  col <- beachcolours(range, ...)
  z <- colourmap(col, range=range)
  return(z)
}

beachcolours <- function(range, sealevel = 0, monochrome=FALSE,
                         ncolours=if(monochrome) 16 else 64,
                         nbeach=1) {
  if(monochrome)
    return(grey(seq(from=0,to=1,length.out=ncolours)))
  stopifnot(is.numeric(range) && length(range) == 2)
  stopifnot(all(is.finite(range)))
  depths <- range[1L]
  peaks <- range[2L]
  dv <- diff(range)/(ncolours - 1L)
  epsilon <- nbeach * dv/2
  lowtide <- max(sealevel - epsilon, depths)
  hightide <-  min(sealevel + epsilon, peaks)
  countbetween <- function(a, b, delta) { max(0, round((b-a)/delta)) }
  nsea <- countbetween(depths, lowtide, dv)
  nbeach <- countbetween(lowtide,  hightide, dv)
  nland <- countbetween(hightide,  peaks, dv)
  colours <- character(0)
  if(nsea > 0)  colours <- rev(rainbow(nsea, start=3/6,end=4/6)) # cyan/blue
  if(nbeach > 0)  colours <- c(colours,
                             rev(rainbow(nbeach, start=3/12,end=5/12))) # green
  if(nland > 0)  colours <- c(colours,
                              rev(rainbow(nland, start=0, end=1/6)))  # red/yellow
  return(colours)
}

