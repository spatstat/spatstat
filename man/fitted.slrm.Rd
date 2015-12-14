\name{fitted.slrm}
\Rdversion{1.1}
\alias{fitted.slrm}
\title{
  Fitted Probabilities for Spatial Logistic Regression
}
\description{
  Given a fitted Spatial Logistic Regression model,
  this function computes the fitted probabilities for each pixel.
}
\usage{
  \method{fitted}{slrm}(object, ...)
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
  This is a method for the generic function
  \code{\link[stats:fitted.values]{fitted}}
  for spatial logistic regression models
  (objects of class \code{"slrm"}, usually obtained
  from the function \code{\link{slrm}}).

  The algorithm computes
  the fitted probabilities of the presence of a random point in each pixel.
}
\value{
  A pixel image (object of class \code{"im"}) containing the
  fitted probability for each pixel.
}
\seealso{
  \code{\link{slrm}},
  \code{\link[stats:fitted.values]{fitted}}
}
\examples{
  X <- rpoispp(42)
  fit <- slrm(X ~ x+y)
  plot(fitted(fit))
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  \email{adrian@maths.uwa.edu.au}
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{models}
\keyword{methods}
