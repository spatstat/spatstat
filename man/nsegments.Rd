\name{nsegments}
\alias{nsegments}
\alias{nsegments.psp}
\title{
   Number of Line Segments in a Line Segment Pattern
}
\description{
   Returns the number of line segments in a line segment pattern.
}
\usage{
   nsegments(x)
   \method{nsegments}{psp}(x)
}
\arguments{
  \item{x}{ A line segment pattern, i.e. an object of class \code{psp}.
}
}
\details{
   This function is generic, but there is at present only one
   method, that for class \code{psp}.
}
\value{
   Integer.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}

\seealso{
  \code{\link{npoints}()}, \code{\link{psp.object}()}
}
\examples{
data(copper)
nsegments(copper$Lines)
nsegments(copper$SouthLines)
}
\keyword{ spatial }
\keyword{ manip }
