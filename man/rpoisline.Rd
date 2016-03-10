\name{rpoisline}
\alias{rpoisline}
\title{Generate Poisson Random Line Process}
\description{
  Generate a random pattern of line segments
  obtained from the Poisson line process.
}
\usage{
 rpoisline(lambda, win=owin())
}
\arguments{
  \item{lambda}{
    Intensity of the Poisson line process.
    A positive number.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
  }
}
\value{
  A line segment pattern (an object of class \code{"psp"}).
}
\details{
  This algorithm generates a realisation
  of the uniform Poisson line process, and clips it to the window
  \code{win}.

  The argument \code{lambda} must be a positive number.
  It controls the intensity of the process. The expected number of
  lines intersecting a convex region of the plane is equal to
  \code{lambda} times the perimeter length of the region.
  The expected total length of the lines crossing a region of the plane
  is equal to \code{lambda * pi} times the area of the region.
}
\seealso{
\code{\link{psp}}
}
\examples{
 # uniform Poisson line process with intensity 10,
 # clipped to the unit square
 rpoisline(10)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{datagen}
