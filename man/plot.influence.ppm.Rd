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

  The display shows circles centred at the data points
  with radii proportional to the influence values.

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
  \adrian
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

