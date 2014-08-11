\name{persp.im}
\alias{persp.im}
\title{Perspective Plot of Pixel Image}
\description{
  Displays a perspective plot of a pixel image.
}
\usage{
  \method{persp}{im}(x, \dots, colmap=NULL)
}
\arguments{
  \item{x}{
    The pixel image to be plotted.
    An object of class \code{"im"} (see \code{\link{im.object}}).
  }
  \item{\dots}{
    Extra arguments passed to \code{\link{persp.default}} to control the
    display.
  }
  \item{colmap}{
    Optional data controlling the colour map. See Details.
  }
}
\value{
  Returns the 3D transformation matrix
  returned by \code{\link{persp.default}}.
}
\details{
  This is the \code{persp} method for the class \code{"im"}.

  The pixel image \code{x} must have real or integer values.
  These values are treated as heights of a surface, and the
  surface is displayed as a perspective plot on the current plot device,
  using equal scales on the \code{x} and \code{y} axes.

  The optional argument \code{colmap} gives an easy way to display
  different altitudes in different colours (if this is what you
  want).
  \itemize{
    \item If \code{colmap} is a colour map (object of class \code{"colourmap"},
    created by the function \code{\link{colourmap}})
    then this colour map will be used to associate
    altitudes with colours.
    \item If \code{colmap} is a character vector, then the range of
    altitudes in the perspective plot will be divided into
    \code{length(colmap)} intervals, and those parts of the surface
    which lie in a particular altitude range will be assigned
    the corresponding colour from \code{colmap}.
    \item If \code{colmap} is a list with entries \code{breaks} and \code{col},
    then \code{colmap$breaks} determines the breakpoints of the altitude
    intervals, and \code{colmap$col} provides the corresponding colours.
  }
  
  Graphical parameters controlling the perspective plot
  are passed through the \code{...} arguments
  directly to the function \code{\link{persp.default}}.
  See the examples in \code{\link{persp.default}} or in
  \code{demo(persp)}.
}
\seealso{
  \code{\link{im.object}},
  \code{\link{plot.im}},
  \code{\link{contour.im}}
}
\examples{
   # an image
   Z <- setcov(owin())
   persp(Z, colmap=terrain.colors(128))
   co <- colourmap(range=c(0,1), col=rainbow(128))
   persp(Z, colmap=co, axes=FALSE, shade=0.3)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
 
 \keyword{hplot}
