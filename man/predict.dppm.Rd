\name{predict.dppm}
\alias{predict.dppm}
\alias{fitted.dppm}
\title{Prediction from a Fitted Determinantal Point Process Model}
\description{
  Given a fitted determinantal point process model,
  these functions compute the fitted intensity.
}
\usage{
  \method{fitted}{dppm}(object, ...)

  \method{predict}{dppm}(object, ...)
}
\arguments{
  \item{object}{
    Fitted determinantal point process model.
    An object of class \code{"dppm"}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{fitted.ppm}} or
    \code{\link{predict.ppm}} respectively.
  }
}
\details{
  These functions are methods for the generic functions
  \code{\link{fitted}} and \code{\link{predict}}.
  The argument \code{object} should be a determinantal point process model
  (object of class \code{"dppm"}) obtained using
  the function \code{\link{dppm}}.

  The \emph{intensity} of the fitted model
  is computed, using \code{\link{fitted.ppm}} or
  \code{\link{predict.ppm}} respectively.
}
\value{
  The value of \code{fitted.dppm} is a numeric vector
  giving the fitted values at the quadrature points.

  The value of \code{predict.dppm} is usually a pixel image
  (object of class \code{"im"}), but see \code{\link{predict.ppm}}
  for details.
}
\seealso{
  \code{\link{dppm}},
  \code{\link{plot.dppm}},
  \code{\link{fitted.ppm}},
  \code{\link{predict.ppm}}
}
\examples{
  fit <- dppm(swedishpines ~ x + y, dppGauss())
  predict(fit)
}
\author{\adrian
  
  ,
  \rolf
  
  and \ege
  
}
\keyword{spatial}
\keyword{models}
