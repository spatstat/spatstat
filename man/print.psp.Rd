\name{print.psp}
\alias{print.psp}
\title{Print Brief Details of a Line Segment Pattern Dataset}
\description{
  Prints a very brief description of a line segment pattern dataset.
}
\usage{
  \method{print}{psp}(x, \dots)
}
\arguments{
  \item{x}{Line segment pattern (object of class \code{"psp"}).}
  \item{\dots}{Ignored.}
}
\details{
  A very brief description of the line segment pattern \code{x} is printed.

  This is a method for the generic function \code{\link{print}}.
}
\seealso{
  \code{\link{print}},
  \code{\link{print.owin}},
  \code{\link{summary.psp}}
}
\examples{
  a <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  a
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{print}

