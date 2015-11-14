\name{print.im}
\alias{print.im}
\title{Print Brief Details of an Image}
\description{
  Prints a very brief description of a pixel image object.
}
\usage{
  \method{print}{im}(x, \dots)
}
\arguments{
  \item{x}{Pixel image (object of class \code{"im"}).}
  \item{\dots}{Ignored.}
}
\details{
  A very brief description of the pixel image \code{x} is printed.

  This is a method for the generic function \code{\link{print}}.
}
\seealso{
  \code{\link{print}},
  \code{\link{im.object}},
  \code{\link{summary.im}}
}
\examples{
  data(letterR)
  U <- as.im(letterR)
  U
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{print}
