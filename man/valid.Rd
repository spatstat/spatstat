\name{valid}
\alias{valid}
\title{
  Check Whether Point Process Model is Valid
}
\description{
  Determines whether a point process model object
  corresponds to a valid point process.  
}
\usage{
  valid(object, \dots)
}
\arguments{
  \item{object}{
    Object of some class, describing a point process model.
  }
  \item{\dots}{
    Additional arguments passed to methods.
  }
}
\details{
  The function \code{valid} is generic,
  with methods for the classes \code{"ppm"} and \code{"dppmodel"}.

  An object representing a point process is called valid if
  all its parameter values are known (for example, no parameter
  takes the value \code{NA} or \code{NaN}) and the parameter values
  correspond to a well-defined point process (for example, the
  parameter values satisfy all the constraints that are imposed by
  mathematical theory.) 
  
  See the methods for further details.
}
\value{
  A logical value, or \code{NA}.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\seealso{
  \code{\link{valid.ppm}},
  \code{\link{valid.detpointprocfamily}}
}
\keyword{spatial}
\keyword{models}
