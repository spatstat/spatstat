\name{rat}
\alias{rat}
\title{
  Ratio object
}
\description{
  Stores the numerator, denominator, and value of a ratio
  as a single object.
}
\usage{
rat(ratio, numerator, denominator, check = TRUE)
}
\arguments{
  \item{ratio,numerator,denominator}{
    Three objects belonging to the same class.
  }
  \item{check}{
    Logical. Whether to check that
    the objects are \code{\link{compatible}}.
  }
}
\details{
  The class \code{"rat"} is a simple mechanism for keeping track of
  the numerator and denominator when calculating a ratio. Its main
  purpose is simply to signal that the object is a ratio.

  The function \code{rat} creates an object of class \code{"rat"}
  given the numerator, the denominator and the ratio.
  No calculation is performed; 
  the three objects are simply stored together.

  The arguments \code{ratio}, \code{numerator}, \code{denominator}
  can be objects of any kind. They should belong to the same class.
  It is assumed that the relationship
  \deqn{
    \mbox{ratio} = \frac{\mbox{numerator}}{\mbox{denominator}}
  }{
    ratio = numerator/denominator
  }
  holds in some version of arithmetic. However, no calculation is
  performed.

  By default the algorithm checks
  whether the three arguments \code{ratio}, \code{numerator},
  \code{denominator} are compatible objects, according to
  \code{\link{compatible}}.
  
  The result is equivalent to \code{ratio} except for the
  addition of extra information.
}
\value{
  An object equivalent to the object \code{ratio}
  except that it also belongs to the class \code{"rat"}
  and has additional attributes \code{numerator} and \code{denominator}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{compatible}}, 
  \code{\link{pool}}
}
\keyword{spatial}
\keyword{manip}
