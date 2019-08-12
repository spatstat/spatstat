\name{data.lppm}
\alias{data.lppm}
\title{Extract Original Data from a Fitted Point Process Model on a Network}
\description{
  Given a fitted point process model on a linear network,
  this function extracts the original point pattern dataset 
  to which the model was fitted.
}
\usage{
  data.lppm(object)
}
\arguments{
  \item{object}{
    fitted point process model
    on a linear network (an object of class \code{"lppm"}).
  }
}
\value{
  A point pattern on a linear network (object of class \code{"lpp"}).
}
\details{
  An object of class \code{"lppm"} represents a point process model
  that has been fitted to a point pattern dataset on a linear network.
  It is typically produced by
  the model-fitting algorithm \code{\link{lppm}}.
  The object contains complete information about the original data
  point pattern to which the model was fitted.
  This function extracts the original data pattern.
}
\seealso{
  \code{\link{lppm}},
  \code{\link{data.ppm}}
}
\examples{
 fit <- lppm(spiders ~ x)
 X <- data.lppm(fit)
 # 'X' is identical to 'spiders'
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
\keyword{models}
