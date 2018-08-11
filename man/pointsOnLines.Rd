\name{pointsOnLines}
\alias{pointsOnLines}
\title{Place Points Evenly Along Specified Lines}
\description{
  Given a line segment pattern, place a series of points
  at equal distances along each line segment.
}
\usage{
pointsOnLines(X, eps = NULL, np = 1000, shortok=TRUE)
}
\arguments{
  \item{X}{A line segment pattern (object of class \code{"psp"}).}
  \item{eps}{Spacing between successive points.}
  \item{np}{Approximate total number of points (incompatible with
    \code{eps}).}
  \item{shortok}{
    Logical. If \code{FALSE}, very short segments
    (of length shorter than \code{eps}) will not generate any points.
    If \code{TRUE}, a very short segment will be represented by its
    midpoint.
  }
}
\details{
  For each line segment in the pattern \code{X}, a succession of points
  is placed along the line segment. These points are equally spaced at
  a distance \code{eps}, except for the first and last points
  in the sequence. 

  The spacing \code{eps} is measured in coordinate units of \code{X}.
  
  If \code{eps} is not given, then it is determined by
  \code{eps = len/np} where \code{len} is the total length of the
  segments in \code{X}. The actual number of points will then be slightly
  larger than \code{np}.
}
\value{
  A point pattern (object of class \code{"ppp"}) in the same window
  as \code{X}. The result also has an attribute called \code{"map"} which
  maps the points to their parent line segments.
}
\seealso{
  \code{\link{psp}},
  \code{\link{ppp}},
  \code{\link{runifpointOnLines}}
}
\examples{
  X <- psp(runif(20), runif(20), runif(20), runif(20),  window=owin())
  Y <- pointsOnLines(X, eps=0.05)
  plot(X, main="")
  plot(Y, add=TRUE, pch="+")
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
