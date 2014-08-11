\name{clip.infline}
\alias{clip.infline}
\title{Intersect Infinite Straight Lines with a Window}
\description{
  Take the intersection between a set of infinite straight lines
  and a window, yielding a set of line segments.
}
\usage{
clip.infline(L, win)
}
\arguments{
  \item{L}{
    Object of class \code{"infline"} specifying a set of infinite
    straight lines in the plane.
  }
  \item{win}{
    Window (object of class \code{"owin"}).
  }
}
\details{
  This function computes the intersection between
  a set of infinite straight lines in the plane
  (stored in an object \code{L} of class \code{"infline"} created by the
  function \code{\link{infline}}) and a window \code{win}.
  The result is a pattern of line segments. 
}
\value{
  A line segment pattern (object of class \code{"psp"}).
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{infline}},\code{\link{psp}}.

  To divide a window into pieces using infinite lines,
  use \code{\link{chop.tess}}.
}
\examples{
  L <- infline(p=1:3, theta=pi/4)
  W <- square(4)
  clip.infline(L, W)
}
\keyword{spatial}
\keyword{math}
