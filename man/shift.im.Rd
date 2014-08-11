\name{shift.im}
\alias{shift.im}
\title{Apply Vector Translation To Pixel Image}
\description{
  Applies a vector shift to a pixel image
}
\usage{
 \method{shift}{im}(X, vec=c(0,0), \dots, origin=NULL)
}
\arguments{
  \item{X}{Pixel image (object of class \code{"im"}).}
  \item{vec}{Vector of length 2 representing a translation.}
  \item{\dots}{Ignored}
  \item{origin}{Character string determining a location
    that will be shifted to the origin. Options are
    \code{"centroid"}, \code{"midpoint"} and \code{"bottomleft"}.
    Partially matched.
  }
}
\value{
  Another pixel image (of class \code{"im"}) representing the
  result of applying the vector shift.
}
\details{
  The spatial location of each pixel in the image
  is translated by the vector \code{vec}.
  This is a method for the generic function \code{\link{shift}}.

  If \code{origin} is given, then it should be one of the character
  strings \code{"centroid"}, \code{"midpoint"} or \code{"bottomleft"}.
  The argument \code{vec} will be ignored; instead the shift will be performed
  so that the specified geometric location is shifted to the origin.
  If \code{origin="centroid"} then the centroid of the image window will be
  shifted to the origin. If \code{origin="midpoint"} then the centre of
  the bounding rectangle of the image will be shifted to the origin.
  If \code{origin="bottomleft"} then the bottom left corner of the
  bounding rectangle of the image will be shifted to the origin.
}
\seealso{
  \code{\link{shift}}
}
\examples{
 # make up an image
 X <- setcov(unit.square())
 plot(X)

 Y <- shift(X, c(10,10))
 plot(Y)
 # no discernible difference except coordinates are different

 shift(X, origin="mid")
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
