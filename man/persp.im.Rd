\name{persp.im}
\alias{persp.im}
\title{Perspective Plot of Pixel Image}
\description{
  Displays a perspective plot of a pixel image.
}
\usage{
  \method{persp}{im}(x, \dots,
                     colmap=NULL, colin=x, apron=FALSE, visible=FALSE)
}
\arguments{
  \item{x}{
    The pixel image to be plotted as a surface.
    An object of class \code{"im"} (see \code{\link{im.object}}).
  }
  \item{\dots}{
    Extra arguments passed to \code{\link{persp.default}} to control the
    display.
  }
  \item{colmap}{
    Optional data controlling the colour map. See Details.
  }
  \item{colin}{
    Optional. Colour input. Another pixel image (of the same dimensions
    as \code{x}) containing the values that will be mapped to colours. 
  }
  \item{apron}{
    Logical. If \code{TRUE}, a grey apron is placed
    around the sides of the perspective plot.
  }
  \item{visible}{
    Logical value indicating whether to compute which pixels of \code{x}
    are visible in the perspective view. See Details.
  }
}
\value{
  (invisibly) the 3D transformation matrix
  returned by \code{\link{persp.default}},
  together with an attribute \code{"expand"} which gives the
  relative scale of the \eqn{z} coordinate.

  If argument \code{visible=TRUE} was given, the return value
  also has an attribute \code{"visible"} which is a pixel image,
  compatible with \code{x}, with logical values which are \emph{TRUE}
  when the corresponding pixel is visible in the perspective view,
  and \code{FALSE} when it is obscured.
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
    \item If \code{colmap} is a function in the \R language
    of the form \code{function(n, ...)}, this function will be called
    with an appropriate value of \code{n} to generate a character
    vector of \code{n} colours. 
      Examples of such functions are
      \code{\link[grDevices]{heat.colors}},
      \code{\link[grDevices]{terrain.colors}},
      \code{\link[grDevices]{topo.colors}} and
      \code{\link[grDevices]{cm.colors}}. 
    \item If \code{colmap} is a function in the \R language of the
    form \code{function(range, ...)} then it will be called
    with \code{range} equal to the range of altitudes,
    to determine the colour values or colour map.
    Examples of such functions are \code{\link{beachcolours}}
    and \code{\link{beachcolourmap}}.
    \item If \code{colmap} is a list with entries \code{breaks} and \code{col},
    then \code{colmap$breaks} determines the breakpoints of the altitude
    intervals, and \code{colmap$col} provides the corresponding colours.
  }
  Alternatively, if the argument \code{colin} (\emph{colour input}) is present,
  then the colour map \code{colmap} will be applied to
  the pixel values of \code{colin} instead of the pixel values of
  \code{x}. The result is a perspective view of a surface with
  heights determined by \code{x} and colours determined by \code{colin}
  (mapped by \code{colmap}).  
  
  If \code{apron=TRUE}, vertical surface is drawn around the boundary
  of the perspective plot, so that the terrain appears to have been
  cut out of a solid material. If colour data were supplied, then
  the apron is coloured light grey.

  Graphical parameters controlling the perspective plot
  are passed through the \code{...} arguments
  directly to the function \code{\link{persp.default}}.
  See the examples in \code{\link{persp.default}} or in
  \code{demo(persp)}.

  The vertical scale is controlled by the argument \code{expand}:
  setting \code{expand=1} will interpret the pixel values as being in
  the same units as the spatial coordinates \eqn{x} and \eqn{y} and
  represent them at the same scale.

  If \code{visible=TRUE}, the algorithm also computes whether each pixel
  in \code{x} is visible in the perspective view. In order to be
  visible, a pixel must not be obscured by another pixel which lies
  in front of it (as seen from the viewing direction), and the
  three-dimensional vector normal to the surface must be pointing toward the
  viewer. The return value of \code{persp.im} then has an attribute
  \code{"visible"} which is a pixel image, compatible with \code{x},
  with pixel value equal to \code{TRUE} if the corresponding pixel in
  \code{x} is visible, and \code{FALSE} if it is not visible.
}
\seealso{
  \code{\link{perspPoints}},
  \code{\link{perspLines}}
  for drawing additional points or lines on the surface.

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

   ## Terrain elevation
   persp(bei.extra$elev, colmap=terrain.colors(128),
         apron=TRUE, theta=-30, phi=20,
         zlab="Elevation", main="", ticktype="detailed",
         expand=6)
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\keyword{spatial}
 
 \keyword{hplot}
