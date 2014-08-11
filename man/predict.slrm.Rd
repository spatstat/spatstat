\name{predict.slrm}
\Rdversion{1.1}
\alias{predict.slrm}
\title{
  Predicted or Fitted Values from Spatial Logistic Regression
}
\description{
  Given a fitted Spatial Logistic Regression model,
  this function computes the fitted probabilities for each pixel,
  or the fitted point process intensity, or the values of the
  linear predictor in each pixel.
}
\usage{
 \method{predict}{slrm}(object, ..., type = "intensity",
                        newdata=NULL, window=NULL)
}
\arguments{
  \item{object}{
    a fitted spatial logistic regression model.
    An object of class \code{"slrm"}.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link[spatstat]{pixellate}}
    determining the pixel resolution for the discretisation
    of the point pattern.
  }
  \item{type}{
    Character string (partially) matching one of 
    \code{"probabilities"}, \code{"intensity"} or \code{"link"}.
  }
  \item{newdata}{
    Optional.  List containing new covariate values for the prediction.
    See Details.
  }
  \item{window}{
    Optional.  New window in which to predict.
    An object of class \code{"owin"}.
  }
}
\details{
  This is a method for \code{\link[stats]{predict}} for spatial logistic
  regression models (objects of class \code{"slrm"}, usually obtained
  from the function \code{\link{slrm}}).

  The argument \code{type} determines which quantity is computed.
  If \code{type="intensity"}), the value of the point process intensity
  is computed at each pixel. If \code{type="probabilities"}) the 
  probability of the presence of a random point in each pixel is
  computed. If \code{type="link"}, the value of the linear predictor is
  computed at each pixel. 

  If \code{newdata = NULL} (the default), the algorithm computes
  fitted values of the model (based on the data that was originally used
  to fit the model \code{object}). 

  If \code{newdata} is given, the algorithm computes predicted values
  of the model, using the new values of the covariates provided by
  \code{newdata}. The argument \code{newdata} should be a list;
  names of entries in the list should correspond
  to variables appearing in the model formula of the \code{object}.
  Each list entry may be a pixel image or a single numeric
  value. 
}
\value{
  A pixel image (object of class \code{"im"}) containing the predicted
  values for each pixel.
}
\seealso{
  \code{\link{slrm}}
}
\examples{
  X <- rpoispp(42)
  fit <- slrm(X ~ x+y)
  plot(predict(fit))

  data(copper)
  X <- copper$SouthPoints
  Y <- copper$SouthLines
  Z <- distmap(Y)
  fitc <- slrm(X ~ Z)
  pc <- predict(fitc)

  Znew <- distmap(copper$Lines)[copper$SouthWindow]
  pcnew <- predict(fitc, newdata=list(Z=Znew))
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
