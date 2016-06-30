\name{plot.laslett}
\alias{plot.laslett}
\title{
  Plot Laslett Transform
}
\description{
  Plot the result of Laslett's Transform.
}
\usage{
\method{plot}{laslett}(x, \dots,
          Xpars = list(box = TRUE, col = "grey"),
          pointpars = list(pch = 3, cols = "blue"),
          rectpars = list(lty = 3, border = "green"))
}
\arguments{
  \item{x}{
    Object of class \code{"laslett"} produced by \code{\link{laslett}}
    representing the result of Laslett's transform.
  }
  \item{\dots}{
    Additional plot arguments passed to \code{\link{plot.solist}}.
  }
  \item{Xpars}{
    A list of plot arguments passed to \code{\link{plot.owin}}
    or \code{\link{plot.im}} to display the original region \code{X}
    before transformation.
  }
  \item{pointpars}{
    A list of plot arguments passed to \code{\link{plot.ppp}}
    to display the tangent points.
  }
  \item{rectpars}{
    A list of plot arguments passed to \code{\link{plot.owin}}
    to display the maximal rectangle.
  }
}
\details{
  This is the \code{plot} method for the class \code{"laslett"}.
  
  The function \code{\link{laslett}} applies Laslett's Transform
  to a spatial region \code{X} and returns an object of class
  \code{"laslett"} representing the result of the transformation.
  The result is plotted by this method.

  The plot function \code{\link{plot.solist}} is used to align
  the before-and-after pictures. See \code{\link{plot.solist}} for
  further options to control the plot.
}
\value{
  None.
}
\author{
  Kassel Hingee and \adrian.
}
\seealso{
  \code{\link{laslett}}
}
\examples{
  b <- laslett(heather$coarse, plotit=FALSE)
  plot(b, main="Heather Data")
}
\keyword{spatial}
\keyword{hplot}

