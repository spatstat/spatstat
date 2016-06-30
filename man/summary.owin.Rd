\name{summary.owin}
\alias{summary.owin}
\title{Summary of a Spatial Window}
\description{
  Prints a useful description of a window object.
}
\usage{
  \method{summary}{owin}(object, \dots)
}
\arguments{
  \item{object}{Window (object of class \code{"owin"}).}
  \item{\dots}{Ignored.}
}
\details{
  A useful description of the window \code{object} is printed.

  This is a method for the generic function \code{\link{summary}}.
}
\seealso{
  \code{\link{summary}},
  \code{\link{summary.ppp}},
  \code{\link{print.owin}}
}
\examples{
  summary(owin())  # the unit square

  data(demopat)
  W <- Window(demopat)  # weird polygonal window
  summary(W)           # describes it

  summary(as.mask(W))  # demonstrates current pixel resolution
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{methods}

