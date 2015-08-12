\name{austates}
\alias{austates}
\docType{data}
\title{
  Australian States and Mainland Territories
}
\description{
  The states and large mainland territories of Australia are
  represented as polygonal regions forming a tessellation.
}
\usage{data(austates)}
\format{
  Object of class \code{"tess"}.
}
\details{
  Western Australia, South Australia, Queensland,
  New South Wales, Victoria and Tasmania (which are states of Australia)
  and the Northern Territory (which is a `territory' of Australia)
  are represented as polygonal regions.

  Offshore territories, and smaller mainland territories, are not shown.

  The dataset \code{austates} is a tessellation object (class \code{"tess"})
  whose tiles are the states and territories.

  The coordinates are latitude and
  longitude in degrees, so the space is effectively a Mercator
  projection of the earth.
}
\source{
  Obtained from the \pkg{oz} package and reformatted.
}
\examples{
data(austates); plot(austates)
}
\keyword{datasets}
