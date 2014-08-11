\name{predict.ppm}
\alias{predict.ppm}
\title{Prediction from a Fitted Point Process Model}
\description{
Given a fitted point process model obtained by \code{\link{ppm}},	
evaluate the spatial trend or the conditional intensity of the model
at new locations.
}
\usage{
   \method{predict}{ppm}(object, window, ngrid=NULL, locations=NULL,
   covariates=NULL, type="trend", X=data.ppm(object),
   \dots, new.coef=NULL, check=TRUE, repair=TRUE)
}
\arguments{
  \item{object}{
    A fitted point process model, typically obtained from
    the model-fitting algorithm \code{\link{ppm}}. An object of
    class \code{"ppm"} (see \code{\link{ppm.object}}).
  }
  \item{window}{
    Optional. A window (object of class \code{"owin"})
    delimiting the locations where predictions
    should be computed. Defaults to the window of the
    original data used to fit the model \code{object}.
  }
  \item{ngrid}{
    Optional. Dimensions of a rectangular grid of locations
    inside \code{window} where the predictions should be computed.
    An integer, or an integer vector of length 2,
    specifying the number of grid points in the \eqn{y} and \eqn{x}
    directions. (Incompatible with \code{locations})
  }
  \item{locations}{
    Optional. Data giving the \eqn{x,y} coordinates (and marks, if required)
    of locations at which predictions should be computed.
    Either a point pattern, or a data frame with columns named \code{x} and
    \code{y}, or a binary image mask. (Incompatible with \code{ngrid})
  }
  \item{covariates}{
    Values of external covariates required by the model.
    Either a data frame or a list of images.
    See Details.
  }
  \item{type}{
    Character string.
    Indicates which property of the fitted model should be predicted.
    Options are \code{"trend"} for the spatial trend, 
    \code{"cif"} or \code{"lambda"} for the conditional intensity,
    and \code{"se"} for the standard error of the fitted spatial trend.
  }
  \item{X}{
    Optional. A point pattern (object of class \code{"ppp"})
    to be taken as the data point pattern when calculating the
    conditional intensity. The default is to use the original data
    to which the model was fitted.
  }
  \item{\dots}{
    Ignored.
  }
  \item{new.coef}{
    Numeric vector of parameter values to replace the 
    fitted model parameters \code{coef(object)}.
  }
  \item{check}{
    Logical value indicating whether to check the internal format
    of \code{object}. If there is any possibility that this object
    has been restored from a dump file, or has otherwise lost track of
    the environment where it was originally computed, set
    \code{check=TRUE}. 
  }
  \item{repair}{
    Logical value indicating whether to repair the internal format
    of \code{object}, if it is found to be damaged. 
  }
}
\value{
  \emph{If \code{locations} is given and is a data frame:}
  a vector of predicted values for the spatial locations
  (and marks, if required) given in \code{locations}.

  \emph{If \code{ngrid} is given, or if \code{locations} is given
    and is a binary image mask:}
  If \code{object} is an unmarked point process,
  the result is a pixel image object (of class \code{"im"}, see
  \code{\link{im.object}}) containing the predictions. 
  If \code{object} is a multitype point process,
  the result is a list of pixel images, containing the predictions
  for each type at the same grid of locations.

  The ``predicted values'' are either values of the spatial trend
  (if \code{type="trend"}), values of the conditional intensity
  (if \code{type="cif"} or \code{type="lambda"}), or
  estimates of standard error for the fitted spatial trend
  (if \code{type="se"}).
}
\details{
  This function computes the spatial trend 
  and the conditional intensity of a fitted spatial point process model,
  and the standard error of the estimate of spatial trend.
  See Baddeley and Turner (2000) for explanation and examples.
 
  Given a point pattern dataset, we may fit
  a point process model to the data using the 
  model-fitting algorithm \code{\link{ppm}}. This
  returns an object of class \code{"ppm"} representing 
  the fitted point process model (see \code{\link{ppm.object}}).
  The parameter estimates in this fitted model can be read off 
  simply by printing the \code{ppm} object.
  The spatial trend and conditional intensity of the 
  fitted model are evaluated using this function \code{predict.ppm}.

  The default action is to create a rectangular grid of points
  in the observation window of the data point pattern, and evaluate
  the spatial trend at these locations.

  The argument \code{type} specifies the values that are computed:
  \describe{
    \item{If \code{type="trend"}:}{
      the ``spatial trend'' of the fitted model is evaluated at each
      required spatial location \eqn{u}.
    }
    \item{If \code{type="cif"}:}{
      the conditional intensity \eqn{\lambda(u, X)}{lambda(u,X)} of the
      fitted model is evaluated at each required spatial location \eqn{u},
      with respect to the data point pattern \eqn{X}.
    }
    \item{If \code{type="se"}:}{
      the estimated (asymptotic) standard error of the
      fitted spatial trend is evaluated at each required spatial location
      \eqn{u}. This is available only for Poisson point process models.
    }
  }
  Note that the ``spatial trend'' is the same as the intensity function
  if the fitted model \code{object} is a Poisson point process.
  However, if the model is not a Poisson process, then the
  ``spatial trend'' is the (exponentiated) first order potential
  and not the intensity of the process. [For example if we fit the
  stationary Strauss process with parameters
  \eqn{\beta}{beta} and \eqn{\gamma}{gamma},
  then the spatial trend is constant and equal to \eqn{\beta}{beta},
  while the intensity is a smaller value that is not easy to
  compute. ]

  The spatial locations where predictions are required,
  are determined by the (incompatible)
  arguments \code{ngrid} and \code{locations}.
  \itemize{
    \item 
    If the argument \code{ngrid} is present, then
    predictions are performed at a rectangular 
    grid of locations in the window \code{window}. 
    The result of prediction will be a pixel image or images.
    \item 
    If \code{locations} is present, then predictions
    will be performed at the spatial locations given by
    this dataset. These may be an arbitrary list of spatial locations,
    or they may be a rectangular grid. 
    The result of prediction will be either a numeric vector
    or a pixel image or images.
    \item 
    If neither \code{ngrid} nor \code{locations} is given, then
    \code{ngrid} is assumed. The value of \code{ngrid} defaults to
    \code{\link{spatstat.options}("npixel")}, which is initialised to 128
    when \pkg{spatstat} is loaded.
  }
  The argument \code{locations} may be a point pattern,
  a data frame or a list specifying arbitrary locations;
  or it may be a binary image mask (an object of class \code{"owin"}
  with type \code{"mask"}) specifying (a subset of) a rectangular
  grid of locations. 

  If \code{locations} is a point pattern (object of class \code{"ppp"}),
  then prediction will be performed at the points of the point pattern.
  The result of prediction will be a vector of predicted values,
  one value for each point.
  If the model is a marked point process, then
  \code{locations} should be a marked point pattern, with marks of the
  same kind as the model; prediction will be performed at these
  marked points.
  The result of prediction will be a vector of predicted values,
  one value for each (marked) point.

  If \code{locations} is a data frame or list, then it must contain
  vectors \code{locations$x} and \code{locations$y} specifying the
  \eqn{x,y} coordinates of the prediction locations. Additionally, if
  the model is a marked point process, then \code{locations} must also contain
  a factor \code{locations$marks} specifying the marks of the
  prediction locations. These vectors must have equal length.
  The result of prediction will be a vector of predicted values,
  of the same length.
  
  If \code{locations} is a binary image mask, then prediction will be
  performed at each pixel in this binary image where the pixel value
  is \code{TRUE} (in other words, at each pixel that is inside the
  window). If the fitted model is an unmarked point process, then the
  result of prediction will be an image. If the fitted model is a
  marked point process, then prediction will
  be performed for each possible value of the mark at each such
  location, and the result of prediction will be a 
  list of images, one for each mark value.
  
  The argument \code{covariates} gives the values of any spatial covariates
  at the prediction locations.
  If the trend formula in the fitted model 
  involves spatial covariates (other than
  the Cartesian coordinates \code{x}, \code{y})
  then \code{covariates} is required.

  The format and use of \code{covariates} are analogous to those of the
  argument of the same name in \code{\link{ppm}}.
  It is either a data frame or a list of images.
  
  If \code{covariates} is a list of images, then
  the names of the entries should correspond to
  the names of covariates in the model formula \code{trend}.
  Each entry in the list must be an image object (of class \code{"im"},
  see \code{\link{im.object}}).
  The software will look up
  the pixel values of each image at the quadrature points.

  If \code{covariates} is a data frame, then the
  \code{i}th row of \code{covariates}
  is assumed to contain covariate data for the \code{i}th location.
  When \code{locations} is a data frame,
  this just means that each row of \code{covariates} contains the
  covariate data for the location specified in the corresponding row of
  \code{locations}. When \code{locations} is a binary image
  mask, the row \code{covariates[i,]} must correspond to the location
  \code{x[i],y[i]} where \code{x = as.vector(raster.x(locations))}
  and \code{y = as.vector(raster.y(locations))}.

  Note that if you only want to use prediction in order to
  generate a plot of the predicted values,
  it may be easier to use \code{\link{plot.ppm}} which calls
  this function and plots the results.
}
\references{
  Baddeley, A. and Turner, R.
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42} (2000) 283--322.
 
  Berman, M. and Turner, T.R. 
  Approximating point process likelihoods with GLIM.
  \emph{Applied Statistics} \bold{41} (1992) 31--38.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{ppm.object}},
  \code{\link{plot.ppm}},
  \code{\link{print.ppm}},
  \code{\link{fitted.ppm}},
  \code{\link{spatstat.options}}
}
\section{Warnings}{
  The current implementation invokes \code{\link{predict.glm}}
  so that \bold{prediction is wrong} if the trend formula in
  \code{object} involves terms in \code{ns()},
  \code{bs()} or \code{poly()}.
  This is a weakness of \code{\link{predict.glm}} itself!
  
  Error messages may be very opaque,
  as they tend to come from deep in the workings of 
  \code{\link{predict.glm}}.
  If you are passing the \code{covariates} argument
  and the function crashes,
  it is advisable to start by checking that all the conditions 
  listed above are satisfied.
}
\examples{
  data(cells)
  \testonly{op <- spatstat.options(npixel=32)}
  m <- ppm(cells, ~ polynom(x,y,2), Strauss(0.05))
  trend <- predict(m, type="trend")
  \dontrun{
  image(trend)
  points(cells)
  }
  cif <- predict(m, type="cif")
  \dontrun{
  persp(cif)
  }
  data(japanesepines)
  mj <- ppm(japanesepines, ~harmonic(x,y,2))
  se <- predict(mj, type="se")

  # prediction at arbitrary locations
  predict(mj, locations=data.frame(x=0.3, y=0.4))

  X <- runifpoint(5, as.owin(japanesepines))
  predict(mj, locations=X)
  predict(mj, locations=X, type="se")

  # multitype
  data(amacrine)
  ma <- ppm(amacrine, ~marks,
     MultiStrauss(c("off","on"),matrix(0.06, 2, 2)))
  Z <- predict(ma)
  Z <- predict(ma, type="cif")
  predict(ma, locations=data.frame(x=0.8, y=0.5,marks="on"), type="cif")

  \testonly{spatstat.options(op)}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
 
 
