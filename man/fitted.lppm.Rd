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
\method{fitted}{lppm}(object, \dots,
                      dataonly = FALSE, new.coef = NULL,
		      leaveoneout = FALSE)
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
  \item{leaveoneout}{
    Logical. If \code{TRUE} the fitted value at each data
    point will be computed using a leave-one-out method. See Details.
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

  If \code{leaveoneout=TRUE}, fitted values will be computed
  for the data points only, using a \sQuote{leave-one-out} rule: 
  the fitted value at \code{X[i]} is effectively computed by
  deleting this point from the data and re-fitting the model to the
  reduced pattern \code{X[-i]}, then predicting the value at
  \code{X[i]}. (Instead of literally performing this calculation,
  we apply a Taylor approximation using the influence function
  computed in \code{\link{dfbetas.ppm}}. 
}
\value{
  A vector containing the values of the fitted spatial trend.

  Entries in this vector correspond to the quadrature points (data or
  dummy points) used to fit the model. 
  The quadrature points can be extracted from \code{object}
  by \code{union.quad(quad.ppm(object))}.
}
\author{
  \adrian
  
  
  \rolf
  
  and \ege
  
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
