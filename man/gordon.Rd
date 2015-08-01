\name{gordon}
\alias{gordon}
\docType{data}
\title{
  People in Gordon Square
}
\description{
  This dataset records the location of
  people sitting on a grass patch in Gordon Square, London,
  at 3pm on a sunny afternoon.

  The dataset \code{gordon} is a point pattern
  (object of class \code{"ppp"}) containing the spatial coordinates
  of each person.

  The grass patch is an irregular polygon with two holes.
  
  Coordinates are given in metres.
}
\usage{data(gordon)}
\examples{
data(gordon)
plot(gordon)
}
\source{
  Andrew Bevan, University College London.
}
\references{
  Baddeley, A., Turner, R., Mateu, J. and Bevan, A. (2013)
  Hybrids of Gibbs point process models and their implementation.
  \emph{Journal of Statistical Software} \bold{55}:11, 1--43.
  \url{http://www.jstatsoft.org/v55/i11/}
}
\keyword{datasets}
\keyword{spatial}
