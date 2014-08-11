\name{is.ppm}
\alias{is.ppm}
\title{Test Whether An Object Is A Fitted Point Process Model}
\description{
  Checks whether its argument is a fitted point process model
  (object of class \code{"ppm"}).
}
\usage{
is.ppm(x)
}
\arguments{
  \item{x}{Any object.}
}
\details{
  This function tests whether the object \code{x}
  is a fitted point process model object of class
  \code{"ppm"}. See \code{\link{ppm.object}} for details
  of this class.
}
\value{
  \code{TRUE} if \code{x} has \code{"ppm"} amongst its classes,
  otherwise \code{FALSE}.
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
\keyword{models}
