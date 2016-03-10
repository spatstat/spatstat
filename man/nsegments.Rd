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
  \item{x}{
    A line segment pattern, i.e. an object of class \code{psp},
    or an object containing a linear network.
  }
}
\details{
   This function is generic, with methods for 
   classes \code{psp}, \code{linnet} and \code{lpp}.
}
\value{
   Integer.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}
}
\seealso{
  \code{\link{npoints}()}, \code{\link{psp.object}()}
}
\examples{
nsegments(copper$Lines)
nsegments(copper$SouthLines)
}
\keyword{spatial}
\keyword{manip}
