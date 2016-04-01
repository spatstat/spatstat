\name{plot.kppm}
\alias{plot.kppm}
\title{Plot a fitted cluster point process}
\description{
  Plots a fitted cluster point process model,
  displaying the fitted intensity and the fitted \eqn{K}-function.
}
\usage{
  \method{plot}{kppm}(x, \dots,
                      what=c("intensity", "statistic", "cluster"),
                      pause=interactive(),
                      xname)
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
  \item{pause}{
    Logical value specifying whether to pause between plots.
  }
  \item{xname}{
    Optional. Character string. The name of the object \code{x}
    for use in the title of the plot.
  }
}
\details{
  This is a method for the generic function
  \code{\link[graphics]{plot}} for the class \code{"kppm"} of fitted
  cluster point process models.
  
  The argument \code{x} should be a cluster point process model
  (object of class \code{"kppm"}) obtained using
  the function \code{\link{kppm}}.

  The choice of plots (and the order in which they are
  displayed) is controlled by the argument \code{what}.
  The options (partially matched) are \code{"intensity"},
  \code{"statistic"} and \code{"cluster"}.

  This command is capable of producing three different plots:
  \describe{

    \item{what="intensity"}{specifies the fitted intensity of the model,
    which is plotted using \code{\link{plot.ppm}}. By default this plot
    is not produced for stationary models.}

    \item{what="statistic"}{specifies the empirical and fitted summary
    statistics, which are plotted using \code{\link{plot.fv}}. This is
    only meaningful if the model has been fitted using the Method of
    Minimum Contrast, and it is turned off otherwise.}

    \item{what="cluster"}{specifies a fitted cluster, which is computed
    by \code{\link{clusterfield}} and plotted by
    \code{\link{plot.im}}. It is only meaningful for Poisson cluster
    (incl. Neyman-Scott) processes, and it is turned off for
    log-Gaussian Cox processes (LGCP). If the model is stationary (and
    non-LGCP) this option is turned on by default and shows a fitted
    cluster positioned at the centroid of the observation window. For
    non-stationary (and non-LGCP) models this option is only invoked if
    explicitly told so, and in that case an additional argument
    \code{locations} (see \code{\link{clusterfield}}) must be given to
    specify where to position the parent point(s) .}
}
Alternatively \code{what="all"} selects all available options.
}
\value{
  Null.
}
\examples{
  data(redwood)
  fit <- kppm(redwood~1, "Thomas")
  plot(fit)
}
\seealso{
  \code{\link{kppm}},
  \code{\link{plot.ppm}},
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\keyword{spatial}
\keyword{models}
