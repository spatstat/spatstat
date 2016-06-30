\name{owin.object}
\alias{owin.object} %DoNotExport
\title{Class owin}
\description{
  A class \code{owin} to define the ``observation window'' of a point pattern
}
\details{
  In the \pkg{spatstat} library, a point pattern dataset must include
  information about the window or region in which the pattern was
  observed. A window is described by an object of class \code{"owin"}.
  Windows of arbitrary shape are supported.
  
  An object of class \code{"owin"} has one of three types:
  \tabular{ll}{
    \code{"rectangle"}: \tab
      a rectangle in the two-dimensional plane with edges parallel to the axes
    \cr
    \code{"polygonal"}: \tab
    a region whose boundary is a polygon or several polygons.
    The region may have holes and may consist of several disconnected pieces.
    \cr
    \code{"mask"}: \tab
      a binary image (a logical matrix)
      set to \code{TRUE} for pixels inside the window and
      \code{FALSE} outside the window.
  }
  Objects of class \code{"owin"} may be created by the function
  \code{\link{owin}}
  and converted from other types of data by the function
  \code{\link{as.owin}}.

  They may be manipulated by the functions 
  \code{\link{as.rectangle}},
  \code{\link{as.mask}}, 
  \code{\link{complement.owin}},
  \code{\link{rotate}},
  \code{\link{shift}},
  \code{\link{affine}},
  \code{\link{erosion}},
  \code{\link{dilation}},
  \code{\link{opening}}
  and
  \code{\link{closing}}.

  Geometrical calculations available for windows include
  \code{\link{area.owin}},
  \code{\link{perimeter}},
  \code{\link{diameter.owin}},
  \code{\link{boundingbox}},
  \code{\link{eroded.areas}},
  \code{\link{bdist.points}},
  \code{\link{bdist.pixels}},
  and
  \code{even.breaks.owin}.
  The mapping between continuous coordinates and pixel raster indices
  is facilitated by the functions
  \code{\link{raster.x}},
  \code{\link{raster.y}} and
  \code{\link{nearest.raster.point}}.

  There is a \code{plot} method for window objects,
  \code{\link{plot.owin}}. This may be useful if you wish to
  plot a point pattern's window without the points for graphical
  purposes.

  There are also methods for
  \code{summary} and \code{print}. 
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{as.rectangle}},
  \code{\link{as.mask}},
  \code{\link{summary.owin}},
  \code{\link{print.owin}},
  \code{\link{complement.owin}},
  \code{\link{erosion}},
  \code{\link{dilation}},
  \code{\link{opening}},
  \code{\link{closing}},
  \code{\link{affine.owin}},
  \code{\link{shift.owin}},
  \code{\link{rotate.owin}},
  \code{\link{raster.x}},
  \code{\link{raster.y}},
  \code{\link{nearest.raster.point}},
  \code{\link{plot.owin}},
  \code{\link{area.owin}},
  \code{\link{boundingbox}},
  \code{\link{diameter}},
  \code{\link{eroded.areas}},
  \code{\link{bdist.points}},
  \code{\link{bdist.pixels}}
}
\section{Warnings}{
  In a window of type \code{"mask"}, the 
  row index corresponds to increasing \eqn{y} coordinate,
  and the column index corresponds to increasing \eqn{x} coordinate.
}
\examples{
 w <- owin()
 w <- owin(c(0,1), c(0,1))
 # the unit square
  
 w <- owin(c(0,1), c(0,2))
 \dontrun{
 if(FALSE) {
   plot(w)
   # plots edges of a box 1 unit x 2 units
   v <- locator() 
   # click on points in the plot window
   # to be the vertices of a polygon 
   # traversed in anticlockwise order 
   u <- owin(c(0,1), c(0,2), poly=v)
   plot(u)
   # plots polygonal boundary using polygon()
   plot(as.mask(u, eps=0.02))
   # plots discrete pixel approximation to polygon
 }
}
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{attribute}
 
 
