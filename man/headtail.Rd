\name{headtail}
\alias{head.ppp}
\alias{head.ppx}
\alias{head.psp}
\alias{head.tess}
\alias{tail.ppp}
\alias{tail.ppx}
\alias{tail.psp}
\alias{tail.tess}
\title{
  First or Last Part of a Spatial Pattern
}
\description{
  Returns the first few elements (\code{head}) or the last few
  elements (\code{tail}) of a spatial pattern.
}
\usage{
  \method{head}{ppp}(x, n = 6L, \dots)

  \method{head}{ppx}(x, n = 6L, \dots)

  \method{head}{psp}(x, n = 6L, \dots)

  \method{head}{tess}(x, n = 6L, \dots)

  \method{tail}{ppp}(x, n = 6L, \dots)

  \method{tail}{ppx}(x, n = 6L, \dots)

  \method{tail}{psp}(x, n = 6L, \dots)

  \method{tail}{tess}(x, n = 6L, \dots)
}
\arguments{
  \item{x}{
    A spatial pattern of geometrical figures,
    such as a spatial pattern of points
    (an object of class \code{"ppp"}, \code{"pp3"}, \code{"ppx"} or
    \code{"lpp"}) or a spatial pattern of line segments
    (an object of class \code{"psp"}) or a tessellation
    (object of class \code{"tess"}).
  }
  \item{n}{
    Integer. The number of elements of the pattern that should be extracted.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  These are methods for the generic functions \code{\link[utils]{head}}
  and \code{\link[utils]{tail}}. They extract the first or last
  \code{n} elements from \code{x} and return them as an object of the
  same kind as \code{x}.

  To inspect the spatial coordinates themselves, use
  \code{\link[utils]{View}(x)}
  or \code{head(as.data.frame(x))}.
}
\value{
  An object of the same class as \code{x}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[utils]{View}}, \code{\link[utils]{edit}}.
  
  Conversion to data frame:
  \code{\link{as.data.frame.ppp}},
  \code{\link{as.data.frame.ppx}},
  \code{\link{as.data.frame.psp}}
}
\examples{
  head(cells)
  tail(as.psp(spiders), 10)
  head(dirichlet(cells), 4)
}
\keyword{spatial}
\keyword{manip}
