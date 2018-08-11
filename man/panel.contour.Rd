\name{panel.contour}
\alias{panel.contour}
\alias{panel.image}
\alias{panel.histogram}
\title{
  Panel Plots using Colour Image or Contour Lines
}
\description{
  These functions can be passed to \code{\link[graphics]{pairs}} or
  \code{\link[graphics]{coplot}}
  to determine what kind of plotting is done in each panel
  of a multi-panel graphical display. 
}
\usage{
panel.contour(x, y, ..., sigma = NULL)

panel.image(x, y, ..., sigma = NULL)

panel.histogram(x, ...)
}
\arguments{
  \item{x,y}{
    Coordinates of points in a scatterplot.
  }
  \item{\dots}{
    Extra graphics arguments, passed to \code{\link{contour.im}},
    \code{\link{plot.im}} or \code{\link[graphics]{rect}}, respectively,
    to control the appearance of the panel.
  }
  \item{sigma}{
    Bandwidth of kernel smoother, on a scale where
    \eqn{x} and \eqn{y} range between 0 and 1.
  }
}
\details{
  These functions can serve as one of the arguments \code{panel},
  \code{lower.panel}, \code{upper.panel}, \code{diag.panel}
  passed to graphics commands like
  \code{\link[graphics]{pairs}} or \code{\link[graphics]{coplot}},
  to determine what kind of plotting is done in each panel
  of a multi-panel graphical display. In particular they work
  with \code{\link{pairs.im}}.

  The functions \code{panel.contour} and \code{panel.contour}
  are suitable for the off-diagonal plots which involve
  two datasets \code{x} and \code{y}.
  They first rescale \code{x} and \code{y} to the unit square,
  then apply kernel smoothing with bandwidth \code{sigma}
  using \code{\link{density.ppp}}.
  Then \code{panel.contour} draws a contour plot
  while \code{panel.image} draws a colour image.

  The function \code{panel.histogram} is suitable for the
  diagonal plots which involve a single dataset \code{x}.
  It displays a histogram of the data.
}
\value{
  Null.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{pairs.im}},
  \code{\link{pairs.default}},
  \code{\link{panel.smooth}}
}
\examples{
  pairs(bei.extra,
        panel      = panel.contour,
        diag.panel = panel.histogram)
  with(bei.extra,
        pairs(grad, elev, 
              panel      = panel.image,
              diag.panel = panel.histogram))
  pairs(marks(finpines), panel=panel.contour, diag.panel=panel.histogram)
}
\keyword{spatial}
\keyword{hplot}

