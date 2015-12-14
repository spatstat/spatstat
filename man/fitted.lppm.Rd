\name{fitted.lppm}
\alias{fitted.lppm}
\title{
  Fitted Intensity for Point Process on Linear Network
}
\description{
  Given a point process model fitted to a point pattern on a linear network,
  compute the fitted intensity of the model
  at the points of the pattern,
  or at the points of the quadrature scheme used to fit the model.
}
\usage{
\method{fitted}{lppm}(object, \dots, dataonly = FALSE, new.coef = NULL)
}
\arguments{
  \item{object}{
    Fitted point process model on a linear network
    (object of class \code{"lppm"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{dataonly}{
    Logical value indicating whether to computed fitted intensities
    at the points of the original point pattern dataset
    (\code{dataonly=TRUE})
    or at all the quadrature points of the quadrature scheme
    used to fit the model (\code{dataonly=FALSE}, the default).
  }
  \item{new.coef}{
    Numeric vector of parameter values to replace the 
    fitted model parameters \code{coef(object)}.
  }
}
\details{
  This is a method for the generic function \code{\link[stats]{fitted}}
  for the class \code{"lppm"} of fitted point process models on a linear
  network.
  
  The locations \eqn{u} at which the fitted conditional intensity/trend
  is evaluated, are the points of the
  quadrature scheme used to fit the model in \code{\link{ppm}}.
  They include the data points (the points of the original point pattern
  dataset \code{x}) and other ``dummy'' points 
  in the window of observation.
}
\value{
  A vector containing the values of the fitted spatial trend.

  Entries in this vector correspond to the quadrature points (data or
  dummy points) used to fit the model. 
  The quadrature points can be extracted from \code{object}
  by \code{union.quad(quad.ppm(object))}.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\seealso{
  \code{\link{lppm}},
  \code{\link{predict.lppm}}
}
\examples{
   fit <- lppm(spiders~x+y)
   a <- fitted(fit)
   b <- fitted(fit, dataonly=TRUE)
}
\keyword{spatial}
\keyword{methods}
\keyword{models}
