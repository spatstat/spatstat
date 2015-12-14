\name{is.ppm}
\alias{is.ppm}
\alias{is.lppm}
\alias{is.kppm}
\alias{is.slrm}
\title{Test Whether An Object Is A Fitted Point Process Model}
\description{
  Checks whether its argument is a fitted point process model
  (object of class \code{"ppm"}, \code{"kppm"}, \code{"lppm"}
  or \code{"slrm"}).
}
\usage{
is.ppm(x)
is.kppm(x)
is.lppm(x)
is.slrm(x)
}
\arguments{
  \item{x}{Any object.}
}
\details{
  These functions test whether the object \code{x}
  is a fitted point process model object of the specified class.

  The result of \code{is.ppm(x)} is \code{TRUE} if \code{x}
  has \code{"ppm"} amongst its classes, and 
  otherwise \code{FALSE}. Similarly for the other functions.
}
\value{
  A single logical value.
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{manip}
\keyword{models}
