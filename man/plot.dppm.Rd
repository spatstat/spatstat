\name{plot.dppm}
\alias{plot.dppm}
\title{Plot a fitted determinantal point process}
\description{
  Plots a fitted determinantal point process model,
  displaying the fitted intensity and the fitted summary function.
}
\usage{
  \method{plot}{dppm}(x, ..., what=c("intensity", "statistic"))
}
\arguments{
  \item{x}{
    Fitted determinantal point process model.
    An object of class \code{"dppm"}.
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
  \code{\link{plot}} for the class \code{"dppm"} of fitted
  determinantal point process models.

  The argument \code{x} should be a determinantal point process model
  (object of class \code{"dppm"}) obtained using
  the function \code{\link{dppm}}.

  The choice of plots (and the order in which they are
  displayed) is controlled by the argument \code{what}.
  The options (partially matched) are \code{"intensity"} and
  \code{"statistic"}.

  This command is capable of producing two different plots:
  \describe{

    \item{what="intensity"}{specifies the fitted intensity of the model,
    which is plotted using \code{\link{plot.ppm}}. By default this plot
    is not produced for stationary models.}

    \item{what="statistic"}{specifies the empirical and fitted summary
    statistics, which are plotted using \code{\link{plot.fv}}. This is
    only meaningful if the model has been fitted using the Method of
    Minimum Contrast, and it is turned off otherwise.}

  }
}
\value{
  Null.
}
\examples{
  fit <- dppm(swedishpines ~ x + y, dppGauss())
  plot(fit)
}
\seealso{
  \code{\link{dppm}},
  \code{\link{plot.ppm}},
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
