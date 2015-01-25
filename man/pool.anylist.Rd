\name{pool.anylist}
\alias{pool.anylist}
\title{
  Pool Data from a List of Objects
}
\description{
  Pool the data from the objects in a list.
}
\usage{
\method{pool}{anylist}(x, ...)
}
\arguments{
  \item{x}{
    A list, belonging to the class \code{"anylist"},
    containing objects that can be pooled.
  }
  \item{\dots}{
    Optional additional objects 
    which can be pooled with the elements of \code{x}.
  }
}
\details{
  The function \code{\link{pool}} is generic. Its purpose is to combine
  data from several objects of the same type (typically computed
  from different datasets) into a common, pooled estimate. 

  The function \code{pool.anyist} is the method
  for the class \code{"anylist"}. It is used when the objects to be
  pooled are given in a list \code{x}.

  Each of the elements of the list \code{x}, and each of the
  subsequent arguments \code{\dots} if provided, must be an object of the same
  class. 
}
\value{
  An object of the same class as each of the entries in \code{x}.
}
\seealso{
  \code{\link{anylist}},
  \code{\link{pool}}.
}
\examples{
   Keach <- anylapply(waterstriders, Kest, ratio=TRUE, correction="iso")
   K <- pool(Keach)
}
\author{
Adrian Baddeley
\email{Adrian.Baddeley@uwa.edu.au}
\url{http://www.maths.uwa.edu.au/~adrian/}
Rolf Turner
\email{r.turner@auckland.ac.nz}
and Ege Rubak
\email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{nonparametric}
