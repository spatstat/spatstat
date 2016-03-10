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
  data(cells)
  D <- density(cells)
  zapsmall.im(D)
}
\author{Ege Rubak \email{rubak@math.aau.dk}
  and Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
}
\keyword{spatial}
\keyword{methods}
\keyword{univar}
