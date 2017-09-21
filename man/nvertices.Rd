\name{nvertices}
\alias{nvertices}
\alias{nvertices.owin}
\alias{nvertices.default}
\title{
  Count Number of Vertices
}
\description{
  Count the number of vertices in an object for which vertices are well-defined.
}
\usage{
  nvertices(x, \dots)

  \method{nvertices}{owin}(x, \dots)

  \method{nvertices}{default}(x, \dots)
}
\arguments{
  \item{x}{
    A window (object of class \code{"owin"}),
    or some other object which has vertices.
  }
  \item{\dots}{
    Currently ignored.
  }
}
\details{
  This function counts the number of vertices of \code{x}
  as they would be returned by \code{\link{vertices}(x)}.
  It is more efficient than executing \code{npoints(vertices(x))}.
}
\value{
  A single integer.
}
\author{
  \spatstatAuthors
  and Suman Rakshit.
}
\seealso{
  \code{\link{vertices}}
}
\examples{
  nvertices(square(2))
  nvertices(letterR)
}
\keyword{spatial}
\keyword{manip}
