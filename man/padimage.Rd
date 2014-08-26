\name{padimage}
\alias{padimage}
\title{
  Pad the Border of a Pixel Image
}
\description{
  Fills the border of a pixel image with a given value or values.
}
\usage{
padimage(X, value, n=1)
}
\arguments{
  \item{X}{
    Pixel image (object of class \code{"im"}).
  }
  \item{value}{
    Single value to be placed around the border of \code{X}.    
  }
  \item{n}{
    Width of border, in pixels.
  }
}
\details{
  The image \code{X} will be expanded by a margin of \code{n}
  pixels, with their values set to \code{value}.

  The argument \code{value} should be a single value (a vector of length
  1), normally a value of the same type as the pixel values of \code{X}.
  It may be \code{NA}. Alternatively if \code{X} is a
  factor-valued image, \code{value} can be one of the levels of
  \code{X}.
}
\value{
  Another object of class \code{"im"}, of the same type as \code{X}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\seealso{
  \code{\link{im}}
}
\examples{
 Z <- setcov(owin())
 plot(padimage(Z, 1, 10))
}
\keyword{spatial}
\keyword{manip}

