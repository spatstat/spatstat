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
  Bevan, A. (2012) Manuscript in preparation.
}
\keyword{datasets}
\keyword{spatial}
