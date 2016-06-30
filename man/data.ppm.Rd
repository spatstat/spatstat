\name{data.ppm}
\alias{data.ppm}
\title{Extract Original Data from a Fitted Point Process Model}
\description{
  Given a fitted point process model,
  this function extracts the original point pattern dataset
  to which the model was fitted.
}
\usage{
  data.ppm(object)
}
\arguments{
  \item{object}{
    fitted point process model (an object of class \code{"ppm"}).
  }
}
\value{
  A point pattern (object of class \code{"ppp"}).
}
\details{
  An object of class \code{"ppm"} represents a point process model
  that has been fitted to data. It is typically produced by
  the model-fitting algorithm \code{\link{ppm}}.
  The object contains complete information about the original data
  point pattern to which the model was fitted.
  This function extracts the original data pattern.

  See \code{\link{ppm.object}} for a list of all operations that can be
  performed on objects of class \code{"ppm"}.
}
\seealso{
  \code{\link{ppm.object}},
  \code{\link{ppp.object}}
}
\examples{
 data(cells)
 fit <- ppm(cells, ~1, Strauss(r=0.1))
 X <- data.ppm(fit)
 # 'X' is identical to 'cells'
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{manip}
\keyword{models}
