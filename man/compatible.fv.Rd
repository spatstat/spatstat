\name{compatible.fv}  %DontDeclareMethods
\alias{compatible.fv}
\title{Test Whether Function Objects Are Compatible}
\description{
  Tests whether two or more function objects (class \code{"fv"})
  are compatible.
}
\usage{
  \method{compatible}{fv}(A, B, \dots)
}
\arguments{
  \item{A,B,\dots}{Two or more function value objects (class \code{"fv"}).}
}
\details{
  An object of class \code{"fv"} is essentially a data frame
  containing several different statistical estimates of the same
  function. Such objects are returned by \code{\link{Kest}} and its
  relatives.

  This command tests whether such objects are compatible
  (so that, for example, they could be added or subtracted).
  It is a method for the generic command \code{\link{compatible}}.
  
  The functions are compatible if they have been evaluated at the
  same sequence of values of the argument \code{r}, and if the
  statistical estimates have the same names.
}
\value{
  Logical value: \code{TRUE} if the objects are compatible, and \code{FALSE}
  if they are not.
}
\seealso{
  \code{\link{eval.fv}} 
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{manip}
