\name{perspPoints}
\alias{perspPoints}
\alias{perspSegments}
\alias{perspLines}
\alias{perspContour}
\title{
  Draw Points or Lines on a Surface Viewed in Perspective
}
\description{
  After a surface has been plotted in a perspective view
  using \code{\link{persp.im}}, these functions can be used to
  draw points or lines on the surface.
}
\usage{
perspPoints(x, y=NULL, \dots, Z, M)

perspLines(x, y = NULL, \dots, Z, M)

perspSegments(x0, y0 = NULL, x1 = NULL, y1 = NULL, \dots, Z, M)

perspContour(Z, M, \dots,
             nlevels=10, levels=pretty(range(Z), nlevels))
}
\arguments{
  \item{x,y}{
    Spatial coordinates,
    acceptable to \code{\link[grDevices]{xy.coords}},
    for the points or lines on the horizontal plane.
  }
  \item{Z}{
    Pixel image (object of class \code{"im"})
    specifying the surface heights.
  }
  \item{M}{
    Projection matrix returned from \code{\link{persp.im}}
    when \code{Z} was plotted.
  }
  \item{\dots}{
    Graphical arguments passed to
    \code{\link[graphics]{points}},
    \code{\link[graphics]{lines}} or
    \code{\link[graphics]{segments}}
    to control the drawing.
  }
  \item{x0,y0,x1,y1}{
    Spatial coordinates of the line segments, on the
    horizontal plane. Alternatively \code{x0} can be a
    line segment pattern (object of class \code{"psp"})
    and \code{y0,x1,y1} can be \code{NULL}.
  }
  \item{nlevels}{Number of contour levels}
  \item{levels}{Vector of heights of contours.}
}
\details{
  After a surface has been plotted in a perspective view,
  these functions can be used to draw points or lines on the surface.

  The user should already have called \code{\link{persp.im}}
  in the form \code{M <- persp(Z, visible=TRUE, ...)} 
  to display the perspective view of the surface \code{Z}.

  Only points and lines which are visible from the viewer's standpoint
  will be drawn.  
}
\value{
  Same as the return value from
  \code{\link[graphics]{points}}
  or \code{\link[graphics]{segments}}. 
}
\seealso{
  \code{\link{persp.im}}
}
\examples{
  M <- persp(bei.extra$elev, colmap=terrain.colors(128),
             apron=TRUE, theta=-30, phi=20,
             zlab="Elevation", main="", 
             expand=6, visible=TRUE, shade=0.3)

  perspContour(bei.extra$elev, M=M, col="pink", nlevels=12)
  perspPoints(bei, Z=bei.extra$elev, M=M, pch=16, cex=0.3, col="chartreuse")
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{hplot}
