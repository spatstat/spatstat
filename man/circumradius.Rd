\name{circumradius}
\alias{circumradius}
\alias{circumradius.owin}
\title{
  Circumradius of a Window
}
\description{
  Compute the circumradius of a window or other spatial object.
}
\usage{
circumradius(x, \dots)

\method{circumradius}{owin}(x, \dots)
}
\arguments{
  \item{x}{
    A window (object of class \code{"owin"}), or another spatial object. 
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} to determine the
    pixel resolution for the calculation.
  }
}
\details{
  The circumradius of a spatial region \eqn{W} is the radius of the
  smallest circle that contains \eqn{W}. 
  
  The function \code{circumradius} is generic, with methods
  for objects of class \code{"owin"} and \code{"linnet"}.
}
\value{
  A single numeric value.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au} 
  
}
\seealso{
 \code{\link{circumradius.linnet}}
}
\examples{
  circumradius(letterR)
}
\keyword{spatial}
\keyword{math}
