\name{as.layered.msr}
\alias{as.layered.msr}
\title{Convert Measure To Layered Object}
\description{
  Converts a measure into a layered object.
}
\usage{
 \method{as.layered}{msr}(X)
}
\arguments{
  \item{X}{
    A measure (object of class \code{"msr"}).
  }
}
\value{
  An object of class \code{"layered"} (see \code{\link{layered}}).
}
\details{
  This function converts the object \code{X}
  into an object of class \code{"layered"}.

  It is a method for the generic \code{\link{as.layered}}
  for the class of measures.

  If \code{X} is a vector-valued measure,
  then \code{as.layered(X)} consists of several layers,
  each containing a scalar-valued measure.
}
\seealso{
  \code{\link{as.layered}},
  \code{\link{msr}}.
}
\examples{
   P <- rpoispp(100)
   fit <- ppm(P ~ x+y)
   rs <- residuals(fit, type="score")
   as.layered(rs)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
 
 
