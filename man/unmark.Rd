\name{unmark}
\alias{unmark}
\alias{unmark.ppp}
\alias{unmark.splitppp}
\alias{unmark.psp}
\alias{unmark.ppx}
\title{Remove Marks}
\description{
  Remove the mark information from a spatial dataset. 
}
\usage{
 unmark(X)
 \method{unmark}{ppp}(X)
 \method{unmark}{splitppp}(X)
 \method{unmark}{psp}(X)
 \method{unmark}{ppx}(X)
}
\arguments{
  \item{X}{A point pattern (object of class \code{"ppp"}),
    a split point pattern (object of class \code{"splitppp"}),
    a line segment pattern (object of class \code{"psp"})
    or a multidimensional space-time point pattern
    (object of class \code{"ppx"}).
  }
}
\value{
  An object of the same class as \code{X}
  with any mark information deleted.
}
\details{
  A `mark' is a value attached to each point in a spatial point pattern,
  or attached to each line segment in a line segment pattern, etc.
  
  The function \code{unmark} is a simple way to 
  remove the marks from such a dataset.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{psp.object}}
}
\examples{
  data(lansing)
  hicks <- lansing[lansing$marks == "hickory", ]
  \dontrun{
  plot(hicks)  # still a marked point pattern, but only 1 value of marks
  plot(unmark(hicks)) # unmarked
  }
}
\author{
  \adrian
  and \rolf
}
\keyword{spatial}
\keyword{manip}
