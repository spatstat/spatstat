\name{plot.influence.ppm}
\alias{plot.influence.ppm}
\title{
  Plot Influence Measure
}
\description{
  Plots an influence measure that has been
  computed by \code{\link{influence.ppm}}.
}
\usage{
 \method{plot}{influence.ppm}(x, ..., multiplot=TRUE)
}
\arguments{
  \item{x}{
    Influence measure (object of class \code{"influence.ppm"}) computed by
    \code{\link{influence.ppm}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.ppp}}
    to control the plotting.
  }
  \item{multiplot}{
    Logical value indicating whether it is permissible to
    plot more than one panel. This happens if the original point process
    model is multitype.
  }
}
\details{
  This is the plot method for objects of class \code{"influence.ppm"}.
  These objects are computed by the command \code{\link{influence.ppm}}.

  For a point process model fitted by maximum likelihood or
  maximum pseudolikelihood (the default), influence values are
  associated with the data points.
  The display shows circles centred at the data points
  with radii proportional to the influence values.
  If the original data were a multitype point pattern, then
  if \code{multiplot=TRUE} (the default), 
  there is one such display for each possible type of point,
  while if \code{multiplot=FALSE} there is a single plot
  combining all data points regardless of type.
  
  For a model fitted by logistic composite likelihood
  (\code{method="logi"} in \code{\link{ppm}}) influence values
  are associated with the data points and also with the
  dummy points used to fit the model. The display consist of two
  panels, for the data points and dummy points respectively,
  showing circles with radii proportional to the influence values.
  If the original data were a multitype point pattern, then
  if \code{multiplot=TRUE} (the default), 
  there is one pair of panels for each possible type of point,
  while if \code{multiplot=FALSE} there is a single plot
  combining all data and dummy points regardless of type.

  Use the argument \code{clipwin} to restrict the plot to a subset
  of the full data.
}
\value{
  None.
}
\references{
  Baddeley, A. and Chang, Y.M. and Song, Y. (2013)
  Leverage and influence diagnostics for spatial point process models.
  \emph{Scandinavian Journal of Statistics} \bold{40}, 86--104.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{influence.ppm}}
}
\examples{
   X <- rpoispp(function(x,y) { exp(3+3*x) })
   fit <- ppm(X, ~x+y)
   plot(influence(fit))
}
\keyword{spatial}
\keyword{models}

