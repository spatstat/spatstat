\name{plot.kppm}
\alias{plot.kppm}
\title{Plot a fitted cluster point process}
\description{
  Plots a fitted cluster point process model,
  displaying the fitted intensity and the fitted \eqn{K}-function.
}
\usage{
  \method{plot}{kppm}(x, ..., what=c("intensity", "statistic"))
}
\arguments{
  \item{x}{
    Fitted cluster point process model.
    An object of class \code{"kppm"}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.ppm}}
    and \code{\link{plot.fv}} to control the plot.
  }
  \item{what}{
    Character vector determining what will be plotted.
  }
}
\details{
  This is a method for the generic function
  \code{\link{plot}} for the class \code{"kppm"} of fitted
  cluster point process models.
  
  The argument \code{x} should be a cluster point process model
  (object of class \code{"kppm"}) obtained using
  the function \code{\link{kppm}}.

  By default, this command will first plot the fitted intensity of the model,
  using \code{\link{plot.ppm}}, and then plot the empirical and fitted
  summary statistics, using \code{\link{plot.fv}}.

  The choice of plots (and the order in which they are
  displayed) is controlled by the argument \code{what}.
  The options (partially matched) are \code{"intensity"} and \code{"statistic"}.

  The option \code{what="intensity"} will be ignored if the
  model is stationary.
}
\value{
  Null.
}
\examples{
  data(redwood)
  fit <- kppm(redwood, ~1, "Thomas")
  plot(fit)
}
\seealso{
  \code{\link{kppm}},
  \code{\link{plot.ppm}},
  \code{\link{plot.minconfit}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
