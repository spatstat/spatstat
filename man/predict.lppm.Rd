\name{predict.lppm}
\alias{predict.lppm}
\title{
  Predict Point Process Model on Linear Network
}
\description{
  Given a fitted point process model on a linear network,
  compute the fitted intensity or conditional intensity of the model.
}
\usage{
\method{predict}{lppm}(object, ...,
                       type = "trend", locations = NULL, new.coef=NULL)
}
\arguments{
  \item{object}{
    The fitted model. An object of class \code{"lppm"},
    see \code{\link{lppm}}.
  }
  \item{type}{
    Type of values to be computed. Either \code{"trend"},
    \code{"cif"} or \code{"se"}.
  }
  \item{locations}{
    Optional. Locations at which predictions should be computed.
    Either a data frame with two columns of coordinates,
    or a binary image mask.
  }
  \item{new.coef}{
    Optional. Numeric vector of model coefficients,
    to be used instead of the fitted coefficients
    \code{coef(object)} when calculating
    the prediction.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link{as.mask}}
    to determine the
    pixel resolution (if \code{locations} is missing).
  }
}
\details{
  This function computes the fitted poin process intensity,
  fitted conditional intensity, or standard error of the fitted
  intensity, for a point process model on a linear network.
  It is a method for the generic \code{\link[stats]{predict}}
  for the class \code{"lppm"}.

  The argument \code{object} should be an object of class \code{"lppm"}
  (produced by \code{\link{lppm}}) representing a point process model
  on a linear network.

  Predicted values are computed at the locations given by the
  argument \code{locations}. If this argument is missing,
  then predicted values are computed at a fine grid of points
  on the linear network.

  \itemize{
    \item 
    If \code{locations} is missing or \code{NULL} (the default),
    the return value is a pixel image (object of class \code{"linim"}
    which inherits class \code{"im"})
    corresponding to a discretisation
    of the linear network, with numeric pixel values giving the
    predicted values at each location on the linear network.
    \item
    If \code{locations} is a data frame, the result is a 
    numeric vector of predicted values at the locations specified by
    the data frame.
    \item
    If \code{locations} is a binary mask, the result is a pixel image
    with predicted values computed at the pixels of the mask.
  }
}
\value{
  A pixel image (object of class \code{"linim"} which inherits
  class \code{"im"}) or
  a numeric vector, depending on the argument \code{locations}.
  See Details.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
}
\seealso{
  \code{\link{lpp}},
  \code{\link{linim}}
}
\examples{
  example(lpp)
  fit <- lppm(X, ~x)
  v <- predict(fit, type="trend")
  plot(v)
}
\references{
  Ang, Q.W. (2010)
  \emph{Statistical methodology for events on a network}.
  Master's thesis, School of Mathematics and Statistics, University of
  Western Australia.
  
  Ang, Q.W., Baddeley, A. and Nair, G. (2012)
  Geometrically corrected second-order analysis of 
  events on a linear network, with applications to
  ecology and criminology.
  \emph{Scandinavian Journal of Statistics} \bold{39}, 591--617.

  McSwiggan, G., Nair, M.G. and Baddeley, A. (2012)
  Fitting Poisson point process models to events 
  on a linear network. Manuscript in preparation.
}
\keyword{spatial}
\keyword{models}
