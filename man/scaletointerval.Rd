\name{scaletointerval}
\alias{scaletointerval}
\alias{scaletointerval.default}
\alias{scaletointerval.im}
\title{Rescale Data to Lie Between Specified Limits}
\description{
  Rescales a dataset so that the values range exactly between the
  specified limits.
}
\usage{
  scaletointerval(x, from=0, to=1)
  \method{scaletointerval}{default}(x, from=0, to=1)
  \method{scaletointerval}{im}(x, from=0, to=1)
}
\arguments{
  \item{x}{Data to be rescaled.}
  \item{from,to}{Lower and upper endpoints of the interval
    to which the values of \code{x} should be rescaled.
  }
}
\details{
  These functions rescale a dataset \code{x}
  so that its values range exactly between the limits
  \code{from} and \code{to}.

  The method for pixel images (objects of class \code{"im"})
  applies this scaling to the pixel values of \code{x}.
  
  Rescaling cannot be performed if the values in \code{x} are not
  interpretable as numeric, or if the values in \code{x} are all equal.
}
\value{
  An object of the same type as \code{x}.
}
\seealso{
  \code{\link{scale}}
}
\examples{
  X <- as.im(function(x,y) {x+y+3}, unit.square())
  summary(X)
  Y <- scaletointerval(X)
  summary(Y)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{univar}
