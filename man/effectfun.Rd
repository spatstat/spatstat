\name{effectfun}
\alias{effectfun}
\title{Compute Fitted Effect of a Spatial Covariate in a Point Process Model}
\description{
  Compute the trend or intensity of a fitted point process model
  as a function of one of its covariates.
}
\usage{
  effectfun(model, covname, ..., se.fit=FALSE)
}
\arguments{
  \item{model}{
    A fitted point process model (object of class
    \code{"ppm"}, \code{"kppm"} or \code{"lppm"}).
  }
  \item{covname}{
    The name of the covariate. A character string.
    (Needed only if the model has more than one covariate.)
  }
  \item{\dots}{
    The fixed values of other covariates (in the form
    \code{name=value}) if required.
  }
  \item{se.fit}{
    Logical. If \code{TRUE}, asymptotic standard errors of the estimates
    will be computed, together with a 95\% confidence interval.
  }
}
\details{
  The object \code{model} should be an object of class
  \code{"ppm"}, \code{"kppm"} or \code{"lppm"}
  representing a point process model fitted to point pattern data.

  The model's trend formula should involve a spatial covariate
  named \code{covname}. This could be \code{"x"} or \code{"y"}
  representing one of the Cartesian coordinates.
  More commonly the covariate
  is another, external variable that was supplied when fitting the model.
    
  The command \code{effectfun} computes the fitted trend 
  of the point process \code{model} as a function of the covariate
  named \code{covname}. 
  The return value can be plotted immediately, giving a
  plot of the fitted trend against the value of the covariate.

  If the model also involves covariates other than \code{covname},
  then these covariates will be held fixed. Values for
  these other covariates must be provided as arguments
  to \code{effectfun} in the form \code{name=value}. 

  If \code{se.fit=TRUE}, the algorithm also calculates
  the asymptotic standard error of the fitted trend,
  and a (pointwise) asymptotic 95\% confidence interval for the
  true trend.
  
  This command is just a wrapper for the prediction method
  \code{\link{predict.ppm}}. For more complicated computations
  about the fitted intensity, use \code{\link{predict.ppm}}.
}
\section{Trend and intensity}{
  For a Poisson point process model, the trend is the same as the
  intensity of the point process. For a more general Gibbs model, the trend
  is the first order potential in the model (the first order term in the
  Gibbs representation). In Poisson or Gibbs models fitted by
  \code{\link{ppm}}, the trend is the only part of the model that
  depends on the covariates.
}
\value{
  A data frame containing a column of values of the covariate and a column
  of values of the fitted trend.
  If \code{se.fit=TRUE}, there are 3 additional columns containing the
  standard error and the upper and lower limits of a confidence interval.

  If the covariate named \code{covname} is numeric (rather than a factor
  or logical variable), the return value is
  also of class \code{"fv"} so that it can be plotted immediately.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{predict.ppm}},
  \code{\link{fv.object}}
}
\examples{
  data(copper)
  X <- copper$SouthPoints
  D <- distmap(copper$SouthLines)
  fit <- ppm(X, ~polynom(Z, 5), covariates=list(Z=D))
  \donttest{
     plot(effectfun(fit, "Z"))
  } 
  plot(effectfun(fit, "Z", se.fit=TRUE), shade=c("hi", "lo"))
  fit <- ppm(X, ~x + polynom(Z, 5), covariates=list(Z=D))
  plot(effectfun(fit, "Z", x=20))
  fit <- ppm(X, ~x)
  plot(effectfun(fit, "x"))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
