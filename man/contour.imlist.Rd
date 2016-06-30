\name{contour.imlist}
\alias{contour.imlist}
\alias{contour.listof}
\title{Array of Contour Plots}
\description{
  Generates an array of contour plots.
}
\usage{
  \method{contour}{imlist}(x, \dots)

  \method{contour}{listof}(x, \dots)
}
\arguments{
  \item{x}{
    An object of the class \code{"imlist"}
    representing a list of pixel images.
    Alternatively \code{x} may belong to the outdated class \code{"listof"}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.solist}} to control the
    spatial arrangement of panels, and arguments passed
    to \code{\link{contour.im}} to control the display of each panel.
  }
}
\value{
  Null.
}
\details{
  This is a method for the generic command
  \code{contour} for the class \code{"imlist"}.
  An object of class \code{"imlist"} represents a list of pixel images.

  (The outdated class \code{"listof"} is also handled.)

  Each entry in the list \code{x} will be displayed as a contour plot,
  in an array of panels laid out on the same graphics display,
  using \code{\link{plot.solist}}. Invididual panels are plotted
  by \code{\link{contour.im}}.
}
\seealso{
  \code{\link{plot.solist}},
  \code{\link{contour.im}}
}
\examples{
# Multitype point pattern
 contour(D <- density(split(amacrine)))
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\keyword{spatial}
\keyword{hplot}
