\name{Extract.im}
\alias{[.im}
\title{Extract Subset of Image}
\description{
  Extract a subset or subregion of a pixel image.
}
\usage{
  \method{[}{im}(x, i, j, \dots, drop=TRUE, tight=FALSE,
                                 raster=NULL, rescue=is.owin(i))
}
\arguments{
  \item{x}{
    A two-dimensional pixel image.
    An object of class \code{"im"}.
  }
  \item{i}{
    Object defining the subregion or subset to be extracted.
    Either a spatial window (an object of class \code{"owin"}), or a
    pixel image with logical values, or a point pattern (an object
    of class \code{"ppp"}), or any type of index that applies to a
    matrix, or something that can be converted to a point pattern
    by \code{\link{as.ppp}} (using the window of \code{x}).
  }
  \item{j}{
    An integer or logical vector serving as the column index if
    matrix indexing is being used. Ignored if \code{i} is a spatial object.
  }
  \item{\dots}{Ignored.}
  \item{drop}{
    Logical value. 
    Locations in \code{w} that lie outside the spatial domain of the
    image \code{x} return a pixel value of \code{NA} if
    \code{drop=FALSE}, and are omitted if \code{drop=TRUE}.
  }
  \item{tight}{
    Logical value. If \code{tight=TRUE}, and if the result of the
    subset operation is an image, the image will be trimmed
    to the smallest possible rectangle.
  }
  \item{raster}{
     Optional. An object of class \code{"owin"} or \code{"im"}
     determining a pixel grid.
  }
  \item{rescue}{
   Logical value indicating whether rectangular blocks of data
   should always be returned as pixel images.
  }
}
\value{
  Either a pixel image or a vector of pixel values. See Details.
}
\details{
  This function extracts a subset of the pixel values in a
  pixel image. (To reassign the pixel values, see \code{\link{[<-.im}}).

  The image \code{x} must be an object of class
  \code{"im"} representing a pixel image defined inside a
  rectangle in two-dimensional space (see \code{\link{im.object}}).

  The subset to be extracted is determined by the arguments \code{i,j}
  according to the following rules (which are checked in this order):
  \enumerate{
    \item
    \code{i} is a spatial object such as a window,
    a pixel image with logical values, or a point pattern;  or
    \item
    \code{i,j} are indices for the matrix \code{as.matrix(x)}; or
    \item
    \code{i} can be converted to a point pattern
    by \code{\link{as.ppp}(i, W=Window(x))},
    and \code{i} is not a matrix.
  }
    
  If \code{i} is a spatial window (an object of class \code{"owin"}),
  the values of the image inside this window are extracted
  (after first clipping the window to the spatial domain of the image
  if necessary).

  If \code{i} is a pixel image with logical values,
  it is interpreted as a spatial window (with \code{TRUE} values
  inside the window and \code{FALSE} outside).

  If \code{i} is a point pattern (an object of class
  \code{"ppp"}), then the values of the pixel image at the points of
  this pattern are extracted. This is a simple way to read the
  pixel values at a given spatial location.

  At locations outside the spatial domain of the image, the pixel
  value is undefined, and is taken to be \code{NA}.  The logical
  argument \code{drop} determines whether such \code{NA} values
  will be returned or omitted.  It also influences the format of
  the return value.

  If \code{i} is a point pattern (or something that can be converted
  to a point pattern), then \code{X[i, drop=FALSE]} is a numeric
  vector containing the pixel values at each of the points of
  the pattern.  Its length is equal to the number of points in the
  pattern \code{i}.  It may contain \code{NA}s corresponding to
  points which lie outside the spatial domain of the image \code{x}.
  By contrast, \code{X[i]} or \code{X[i, drop=TRUE]} contains only
  those pixel values which are not \code{NA}. It may be shorter.

  If \code{i} is a spatial window then \code{X[i, drop=FALSE]} is
  another pixel image of the same dimensions as \code{x} obtained
  by setting all pixels outside the window \code{i} to have value
  \code{NA}. When the result is displayed by \code{\link{plot.im}}
  the effect is that the pixel image \code{x} is clipped to the
  window \code{i}.

  If \code{i} is a spatial window then \code{X[i, drop=TRUE]}
  is either:
  \itemize{
    \item a numeric vector containing the pixel values for all pixels
          that lie inside the window \code{i}.
          This happens if \code{i} is \emph{not} a rectangle
          (i.e. \code{i$type != "rectangle"})
          or if \code{rescue=FALSE}.
    \item  a pixel image.
          This happens only if
          \code{i} is a rectangle (\code{i$type = "rectangle"})
          and \code{rescue=TRUE} (the default).
  }

  If the optional argument \code{raster} is given, then it should
  be a binary image mask or a pixel image. Then
  \code{x} will first be converted to an image defined on the
  pixel grid implied by \code{raster}, before the subset operation
  is carried out.
  In particular, \code{x[i, raster=i, drop=FALSE]} will return
  an image defined on the same pixel array as the object \code{i}.

  If \code{i} does not satisfy any of the conditions above, then
  the algorithm attempts to interpret \code{i} and \code{j}
  as indices for the matrix \code{as.matrix(x)}.
  Either \code{i} or \code{j} may be missing or blank.
  The result is usually a vector or matrix of pixel values.
  Exceptionally the result is a pixel image if \code{i,j} determines
  a rectangular subset of the pixel grid, and if the user specifies
  \code{rescue=TRUE}.

  Finally, if none of the above conditions is met,
  the object \code{i} may also be a data frame or list of \code{x,y}
  coordinates which will be converted to a point pattern, taking the
  observation window to be \code{Window(x)}. Then the pixel values
  at these points will be extracted as a vector.
}
\section{Warnings}{
  If you have a 2-column matrix containing the \eqn{x,y} coordinates
  of point locations, then to prevent this being interpreted as an
  array index, you should convert it to a \code{data.frame}
  or to a point pattern.

  If \code{W} is a window or a pixel image, then \code{x[W, drop=FALSE]} 
  will return an image defined on the same pixel array
  as the original image \code{x}. If you want to obtain an image
  whose pixel dimensions agree with those of \code{W}, use the
  \code{raster} argument, \code{x[W, raster=W, drop=FALSE]}.
}
\seealso{
  \code{\link{im.object}},
  \code{\link{[<-.im}},
  \code{\link{ppp.object}},
  \code{\link{as.ppp}},
  \code{\link{owin.object}},
  \code{\link{plot.im}}
}
\examples{
 # make up an image
 X <- setcov(unit.square())
 plot(X)

 # a rectangular subset
 W <- owin(c(0,0.5),c(0.2,0.8))
 Y <- X[W]
 plot(Y)

 # a polygonal subset
 R <- affine(letterR, diag(c(1,1)/2), c(-2,-0.7))
 plot(X[R, drop=FALSE])
 plot(X[R, drop=FALSE, tight=TRUE])

 # a point pattern
 P <- rpoispp(20)
 Y <- X[P]

 # look up a specified location
 X[list(x=0.1,y=0.2)]

 # 10 x 10 pixel array
 X <- as.im(function(x,y) { x + y }, owin(c(-1,1),c(-1,1)), dimyx=10)
 # 100 x 100 
 W <- as.mask(disc(1, c(0,0)), dimyx=100)
 # 10 x 10 raster
 X[W,drop=FALSE]
 # 100 x 100 raster
 X[W, raster=W, drop=FALSE]
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
