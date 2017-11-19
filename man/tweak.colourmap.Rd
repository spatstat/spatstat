\name{tweak.colourmap}
\alias{tweak.colourmap}
\title{
  Change Colour Values in a Colour Map
}
\description{
  Assign new colour values to some of the entries in a colour map.
}
\usage{
tweak.colourmap(m, col, ..., inputs=NULL, range=NULL)
}

\arguments{
  \item{m}{
    A colour map (object of class \code{"colourmap"}).
  }
  \item{inputs}{
    Input values to the colour map, to be assigned new colours.
    Incompatible with \code{range}.
  }
  \item{range}{
    Numeric vector of length 2 specifying a range of numerical values
    which should be assigned a new colour.
    Incompatible with \code{inputs}.
  }
  \item{col}{
    Replacement colours for the specified \code{inputs} or
    the specified \code{range} of values.
  }
  \item{\dots}{Other arguments are ignored.}
}
\details{
  This function changes the colour map \code{m}
  by assigning new colours to each of the
  input values specified by \code{inputs},
  or by assigning a single new colour
  to the range of input values specified by \code{range}.
  
  The modified colour map is returned.
}
\value{
  Another colour map (object of class \code{"colourmap"}).
}
\seealso{
  \code{\link{colourmap}}, 
  \code{\link{interp.colourmap}}, 
  \code{\link{colouroutputs}}, 
  \code{\link[spatstat:colourtools]{colourtools}}.
}
\examples{
  co <- colourmap(rainbow(32), range=c(0,1))
  plot(tweak.colourmap(co, inputs=c(0.5, 0.6), "white"))
  plot(tweak.colourmap(co, range=c(0.5,0.6), "white"))
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{color}

