\name{colouroutputs}
\alias{colouroutputs}
\alias{colouroutputs<-}
\title{
  Extract or Assign Colour Values in a Colour Map
}
\description{
  Extract the colour values in a colour map,
  or assign new colour values.
}
\usage{
colouroutputs(x)

colouroutputs(x) <- value
}

\arguments{
  \item{x}{
    A colour map (object of class \code{"colourmap"}).
  }
  \item{value}{
    A vector of values that can be interpreted as colours.
  }
}
\details{
  An object of class \code{"colourmap"} is effectively a function that
  maps its inputs (numbers or factor levels) to colour values.
  
  The command \code{colouroutputs(x)} extracts the colour values
  in the colour map \code{x}.

  The assignment \code{colouroutputs(x) <- value} replaces the colour
  values in the colour map \code{x} by the entries in 
  \code{value}. The replacement vector \code{value} should have the
  same length as \code{colouroutputs(x)}, and its entries should be
  interpretable as colours.

  To change only some of the colour values in a colour map,
  it may be easier to use \code{\link{tweak.colourmap}}.
}
\value{
  The result of \code{colouroutputs} is a character vector of
  colour values. The result of the assignment
  \code{colouroutputs(x) <- value}  is 
  another colour map (object of class \code{"colourmap"}).
}
\seealso{
  \code{\link{colourmap}}, 
  \code{\link{interp.colourmap}}, 
  \code{\link{tweak.colourmap}}, 
  \code{\link[spatstat:colourtools]{colourtools}}.
}
\examples{
  m <- colourmap(rainbow(5), range=c(0,1))
  m
  # reverse order of colours
  colouroutputs(m) <- rev(colouroutputs(m))
  m
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{color}




