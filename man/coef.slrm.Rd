\name{coef.slrm}
\Rdversion{1.1}
\alias{coef.slrm}
\title{
  Coefficients of Fitted Spatial Logistic Regression Model 
}
\description{
  Extracts the coefficients (parameters) from a fitted
  Spatial Logistic Regression model.
}
\usage{
  \method{coef}{slrm}(object, ...)
}

\arguments{
  \item{object}{
    a fitted spatial logistic regression model.
    An object of class \code{"slrm"}.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a method for \code{\link{coef}} for fitted spatial logistic
  regression models (objects of class \code{"slrm"}, usually obtained
  from the function \code{\link{slrm}}).

  It extracts the fitted canonical parameters, i.e.\ the coefficients in the
  linear predictor of the spatial logistic regression.
}
\value{
  Numeric vector of coefficients.
}
\seealso{
  \code{\link{slrm}}
}
\examples{
  X <- rpoispp(42)
  fit <- slrm(X ~ x+y)
  coef(fit)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
\keyword{methods}
