\name{identify.ppp}
\alias{identify.ppp}
\title{Identify Points in a Point Pattern}
\description{
  If a point pattern is plotted in the graphics window,
  this function will find the point of the pattern which is nearest to
  the mouse position, and print its mark value (or its serial number
  if there is no mark).
}
\usage{
  \method{identify}{ppp}(x, \dots)
}
\arguments{
  \item{x}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[graphics]{identify.default}}.
  }
}
\value{
  If \code{x} is unmarked, the result is 
  a vector containing the serial numbers of the points in the pattern
  \code{x} that were identified.
  If \code{x} is marked, the result is a 
  2-column matrix, the first column containing the serial numbers
  and the second containing the marks for these points.
}
\details{
  This is a method for the generic function \code{\link[graphics]{identify}}
  for point pattern objects.

  The point pattern \code{x} should first be plotted
  using \code{\link{plot.ppp}}. Then \code{identify(x)}
  reads the position of the graphics pointer each time the
  left mouse button is pressed.  It then finds 
  the point of the pattern \code{x} closest to the mouse position.
  If this closest point is sufficiently close to the mouse pointer,
  its index (and its mark if any) 
  will be returned as part of the value of the call.

  Each time a point of the pattern is identified,
  text will be displayed next to the point,
  showing its serial number (if \code{x} is unmarked)
  or its mark value (if \code{x} is marked).
}
\seealso{
  \code{\link[graphics]{identify}},
  \code{\link{clickppp}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{iplot}
