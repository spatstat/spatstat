\name{owin}
\alias{owin}
\title{Create a Window}
\description{
  Creates an object of class \code{"owin"} representing 
  an observation window in the two-dimensional plane
}
\usage{
  owin(xrange=c(0,1), yrange=c(0,1), ..., poly=NULL, mask=NULL,
unitname=NULL, xy=NULL)
}
\arguments{
  \item{xrange}{\eqn{x} coordinate limits of enclosing box}
  \item{yrange}{\eqn{y} coordinate limits of enclosing box}
  \item{\dots}{Ignored.}
  \item{poly}{
    Optional. Polygonal boundary of window.
    Incompatible with \code{mask}.
  }
  \item{mask}{
    Optional. Logical matrix giving binary image of window.
    Incompatible with \code{poly}.
  }
  \item{unitname}{
    Optional. Name of unit of length. Either a single character string,
    or a vector of two character strings giving the
    singular and plural forms, respectively.
  }
  \item{xy}{
    Optional. List with components \code{x} and \code{y} specifying the
    pixel coordinates for \code{mask}.
  }
}
\value{
An object of class \code{"owin"} 
describing a window in the two-dimensional plane.
}
\details{
  In the \pkg{spatstat} library, a point pattern dataset must include
  information about the window of observation. This is represented by
  an object of class \code{"owin"}.
  See \code{\link{owin.object}} for an overview.

  To create a window in its own right, 
  users would normally invoke \code{owin},
  although sometimes \code{\link{as.owin}} may be convenient.

  A window may be rectangular, polygonal, or a mask (a binary image).
  \itemize{
    \item
    \bold{rectangular windows:}
    If only \code{xrange} and \code{yrange} are given, then
    the window will be rectangular, with its \eqn{x} and \eqn{y}
    coordinate dimensions given by these two arguments
    (which must be vectors of length 2).
    If no arguments are given at all, the default is the unit square
    with dimensions \code{xrange=c(0,1)} and \code{yrange=c(0,1)}.
    \item
    \bold{polygonal windows:}
    If \code{poly} is given, then the window will be polygonal.
    \itemize{
      \item
      \emph{single polygon:}
      If \code{poly} is a matrix or data frame with two columns, or
      a structure with two component vectors \code{x} and \code{y}
      of equal length,
      then these values are interpreted as the cartesian coordinates
      of the vertices of a polygon circumscribing the window.
      The vertices must be listed \emph{anticlockwise}.
      No vertex should be repeated (i.e. do not repeat the first
      vertex).
      \item
      \emph{multiple polygons or holes:}
      If \code{poly} is a list, each entry \code{poly[[i]]}
      of which is a matrix or data frame with two columns
      or a structure with
      two component vectors \code{x} and \code{y} of equal length,
      then the successive list members \code{poly[[i]]}
      are interpreted as separate polygons which together
      make up the boundary of the window.
      The vertices of each polygon must be listed \emph{anticlockwise}
      if the polygon is part of the external boundary,
      but \emph{clockwise}
      if the polygon is the boundary of a hole in the window.
      Again, do not repeat any vertex.
    }
    \item
    \bold{binary masks:}
    If \code{mask} is given, then the window will be a binary image.
    \itemize{
      \item
      \emph{Specified by logical matrix:}
      Normally the argument \code{mask} should be a logical matrix
      such that \code{mask[i,j]} is \code{TRUE} if the point
      \code{(x[j],y[i])} belongs to the window, and \code{FALSE} if it
      does not. Note carefully that rows of \code{mask} correspond to the
      \eqn{y} coordinate, and columns to the \eqn{x} coordinate.
      Here \code{x} and \code{y} are vectors of \eqn{x} and \eqn{y}
      coordinates equally spaced over \code{xrange} and \code{yrange}
      respectively. The pixel coordinate vectors \code{x} and \code{y}
      may be specified explicitly using the argument \code{xy}, which
      should be a list containing components \code{x} and \code{y}.
      Alternatively there is a sensible default.
      \item
      \emph{Specified by list of pixel coordinates:}
      Alternatively the argument \code{mask} can be a data frame
      with 2 or 3 columns. If it has 2 columns, it is expected to
      contain the spatial coordinates of all the
      pixels which are inside the window.
      If it has 3 columns,
      it should contain the spatial coordinates \eqn{(x,y)}
      of every pixel in the grid, and the logical value associated
      with each pixel. The pixels may be listed in any order.
    }
  }
  To create a window which is mathematically
  defined by inequalities in the Cartesian coordinates, 
  use \code{\link{raster.x}()} and \code{\link{raster.y}()}
  as in the examples below.

  Functions \code{\link{square}} and \code{\link{disc}}
  will create square and circular windows, respectively.
}
\section{Validity of polygon data}{
  Polygon data may contain geometrical inconsistencies such as
  self-intersections and overlaps. These inconsistencies must be
  removed to prevent problems in other \pkg{spatstat} functions.
  By default, polygon data will be repaired automatically
  using polygon-clipping code.
  The repair process may change the number of vertices in a polygon
  and the number of polygon components. 
  To disable the repair process, set \code{spatstat.options(fixpolygons=FALSE)}.
}
\seealso{
  \code{\link{square}},
  \code{\link{disc}},
  \code{\link{owin.object}},
  \code{\link{as.owin}},
  \code{\link{complement.owin}},
  \code{\link{ppp.object}},
  \code{\link{ppp}}
}
\examples{
  w <- owin()
  w <- owin(c(0,1), c(0,1))
  # the unit square

  w <- owin(c(10,20), c(10,30), unitname=c("foot","feet"))
  # a rectangle of dimensions 10 x 20 feet
  # with lower left corner at (10,10)

  # polygon (diamond shape)
  w <- owin(poly=list(x=c(0.5,1,0.5,0),y=c(0,1,2,1)))
  w <- owin(c(0,1), c(0,2), poly=list(x=c(0.5,1,0.5,0),y=c(0,1,2,1)))

  # polygon with hole
  ho <- owin(poly=list(list(x=c(0,1,1,0), y=c(0,0,1,1)),
                       list(x=c(0.6,0.4,0.4,0.6), y=c(0.2,0.2,0.4,0.4))))
  
  w <- owin(c(-1,1), c(-1,1), mask=matrix(TRUE, 100,100))
          # 100 x 100 image, all TRUE
  X <- raster.x(w)
  Y <- raster.y(w)
  wm <- owin(w$xrange, w$yrange, mask=(X^2 + Y^2 <= 1))
          # discrete approximation to the unit disc

  \dontrun{
  plot(c(0,1),c(0,1),type="n")
  bdry <- locator()
  # click the vertices of a polygon (anticlockwise)
  }
  \testonly{
  bdry <- list(x=c(0.1,0.3,0.7,0.4,0.2),
               y=c(0.1,0.1,0.5,0.7,0.3))
  }
  w <- owin(poly=bdry)
  \dontrun{plot(w)}
 
 \dontrun{
 im <- as.logical(matrix(scan("myfile"), nrow=128, ncol=128))
 # read in an arbitrary 128 x 128 digital image from text file
 rim <- im[, 128:1]
 # Assuming it was given in row-major order in the file
 # i.e. scanning left-to-right in rows from top-to-bottom,
 # the use of matrix() has effectively transposed rows & columns,
 # so to convert it to our format just reverse the column order.
 w <- owin(mask=rim)
 plot(w)
 # display it to check!
 }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
