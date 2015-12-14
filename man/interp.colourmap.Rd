\name{interp.colourmap}
\alias{interp.colourmap}
\title{
  Interpolate smoothly between specified colours
}
\description{
  Given a colourmap object which maps numbers to colours,
  this function interpolates smoothly between the colours,
  yielding a new colour map.
}
\usage{
interp.colourmap(m, n = 512)
}

\arguments{
  \item{m}{
    A colour map (object of class \code{"colourmap"}).
  }
  \item{n}{
    Number of colour steps to be created in the new colour map.
  }
}
\details{
  Given a colourmap object \code{m}, which maps numerical values to
  colours, this function interpolates the mapping, 
  yielding a new colour map.

  This makes it easy to build a colour map that
  has smooth gradation between different colours or shades.
  First specify a small vector of numbers \code{x} which should be mapped to
  specific colours \code{y}. Use \code{m <- colourmap(y, inputs=x)}
  to create a colourmap that represents this simple
  mapping. Then apply \code{interp.colourmap(m)}
  to obtain a smooth transition between these points.
}
\value{
  Another colour map (object of class \code{"colourmap"}).
}
\seealso{
  \code{\link{colourmap}}, 
  \code{\link{tweak.colourmap}}, 
  \code{\link[spatstat:colourtools]{colourtools}}.
}
\examples{
  co <- colourmap(inputs=c(0, 0.5, 1), c("black", "red", "white"))
  plot(interp.colourmap(co))
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{color}

