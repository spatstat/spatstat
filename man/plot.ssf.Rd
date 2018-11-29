\name{plot.ssf}
\alias{plot.ssf}
\alias{image.ssf}
\alias{contour.ssf}
\title{
  Plot a Spatially Sampled Function
}
\description{
  Plot a spatially sampled function object.
}
\usage{
\method{plot}{ssf}(x, \dots,
                   how = c("smoothed", "nearest", "points"),
                   style = c("image", "contour", "imagecontour"),
                   sigma = NULL, contourargs=list())

\method{image}{ssf}(x, \dots)

\method{contour}{ssf}(x, ..., main, sigma = NULL)
}
\arguments{
  \item{x}{
    Spatially sampled function (object of class \code{"ssf"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[graphics]{image.default}}
    or 
    \code{\link[spatstat]{plot.ppp}} to control the plot.
  }
  \item{how}{
    Character string determining whether to display the
    function values at the data points (\code{how="points"}),
    a smoothed interpolation of the function
    (\code{how="smoothed"}), or the function value at the
    nearest data point (\code{how="nearest"}).
  }
  \item{style}{
    Character string indicating whether to plot the smoothed function as
    a colour image, a contour map, or both.
  }
  \item{contourargs}{
    Arguments passed to \code{\link[graphics]{contour.default}}
    to control the contours, if \code{style="contour"} or
    \code{style="imagecontour"}.
  }
  \item{sigma}{
    Smoothing bandwidth for smooth interpolation.
  }
  \item{main}{
    Optional main title for the plot.
  }
}
\details{
  These are methods for the generic
  \code{\link[graphics]{plot}},
  \code{\link[graphics]{image}} and
  \code{\link[graphics]{contour}}
  for the class \code{"ssf"}.
  
  An object of class \code{"ssf"} represents a
  function (real- or vector-valued) that has been
  sampled at a finite set of points.

  For \code{plot.ssf} there are three types of display.
  If \code{how="points"} the exact function values
  will be displayed as circles centred at the locations where they
  were computed. If \code{how="smoothed"} (the default) these
  values will be kernel-smoothed using \code{\link{Smooth.ppp}}
  and displayed as a pixel image.
  If \code{how="nearest"} the values will be interpolated
  by nearest neighbour interpolation using \code{\link{nnmark}}
  and displayed as a pixel image.

  For \code{image.ssf} and \code{contour.ssf} the values are
  kernel-smoothed before being displayed.
}
\value{
  \code{NULL}.
}
\references{
  Baddeley, A. (2017)
  Local composite likelihood for spatial point processes.
  \emph{Spatial Statistics} \bold{22}, 261--295.

  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  Chapman and Hall/CRC Press.
}
\author{
  \adrian.
}
\seealso{
\code{\link{ssf}}
}
\examples{
  a <- ssf(cells, nndist(cells, k=1:3))
  plot(a, how="points")
  plot(a, how="smoothed")
  plot(a, how="nearest")
}
\keyword{spatial}
\keyword{hplot}
