\name{zapsmall.im}
\alias{zapsmall.im}
\title{Rounding of Pixel Values}
\description{
  Modifies a pixel image, identifying those pixels that have values 
  very close to zero, and replacing the value by zero.
}
\usage{
  zapsmall.im(x, digits)
}
\arguments{
  \item{x}{Pixel image (object of class \code{"im"}).}
  \item{digits}{
    Argument passed to \code{\link{zapsmall}}
    indicating the precision to be used.
  }
}
\details{
  The function \code{\link{zapsmall}} is applied to each pixel value
  of the image \code{x}.
}
\value{
  Another pixel image.
}
\seealso{
  \code{\link{zapsmall}}
}
\examples{
  Z <- as.im(function(x,y) { exp(-40*(x+y)) }, square(1), dimyx=32)
  zapsmall.im(Z)
}
\author{\ege
  and \adrian
}
\keyword{spatial}
\keyword{methods}
\keyword{univar}
