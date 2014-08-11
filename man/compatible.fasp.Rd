\name{compatible.fasp}  %DontDeclareMethods
\alias{compatible.fasp}
\title{Test Whether Function Arrays Are Compatible}
\description{
  Tests whether two or more function arrays (class \code{"fasp"})
  are compatible.
}
\usage{
  \method{compatible}{fasp}(A, B, \dots)
}
\arguments{
  \item{A,B,\dots}{Two or more function arrays (object of class \code{"fasp"}).}
}
\details{
  An object of class \code{"fasp"} can be regarded as an array
  of functions. Such objects are returned by the
  command \code{\link{alltypes}}.

  This command tests whether such objects are compatible
  (so that, for example, they could be added or subtracted).
  It is a method for the generic command \code{\link{compatible}}.
  
  The function arrays are compatible if the arrays have the same dimensions,
  and the corresponding elements in each cell of the array
  are compatible as defined by \code{\link{compatible.fv}}.
}
\value{
  Logical value: \code{TRUE} if the objects are compatible, and \code{FALSE}
  if they are not.
}
\seealso{
  \code{\link{eval.fasp}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
