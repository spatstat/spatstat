\name{compatible.fv}  
\alias{compatible.fv}
\title{Test Whether Function Objects Are Compatible}
\description{
  Tests whether two or more function objects (class \code{"fv"})
  are compatible.
}
\usage{
  \method{compatible}{fv}(A, B, \dots, samenames=TRUE)
}
\arguments{
  \item{A,B,\dots}{Two or more function value objects (class
    \code{"fv"}).}
  \item{samenames}{
    Logical value indicating whether to check for complete agreement
    between the column names of the objects (\code{samenames=TRUE}, the
    default) or just to check that the name of the function argument
    is the same (\code{samenames=FALSE}).
  }
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
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{manip}
