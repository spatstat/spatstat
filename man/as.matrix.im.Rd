\name{as.matrix.im}
\alias{as.matrix.im}
\alias{as.array.im}
\title{Convert Pixel Image to Matrix or Array}
\description{
  Converts a pixel image to a matrix or an array.
}
\usage{
  \method{as.matrix}{im}(x, ...)
  \method{as.array}{im}(x, ...)
}
\arguments{
  \item{x}{A pixel image (object of class \code{"im"}).}
  \item{\dots}{See below.}
}
\details{
  The function \code{as.matrix.im} converts the pixel image \code{x}
  into a matrix containing the pixel values.
  It is handy when you want to extract a summary of the pixel values.
  See the Examples.
  
  The function \code{as.array.im} converts the pixel image to an array.
  By default this is a three-dimensional array of dimension \eqn{n} by
  \eqn{m} by \eqn{1}. 
  If the extra arguments \code{\dots} are given, they will be passed
  to \code{\link{array}}, and they may change the dimensions of the
  array.
}
\value{
  A matrix or array.
}
\seealso{
  \code{\link{as.matrix.owin}}
}
\examples{
   # artificial image
   Z <- setcov(square(1))

   M <- as.matrix(Z)

   median(M)
   
   \dontrun{
   # plot the cumulative distribution function of pixel values
   plot(ecdf(as.matrix(Z)))
   }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
