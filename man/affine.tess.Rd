\name{affine.tess} 
\alias{reflect.tess} 
\alias{shift.tess} 
\alias{rotate.tess} 
\alias{scalardilate.tess} 
\alias{affine.tess} 
\title{Apply Geometrical Transformation To Tessellation}
\description{
  Apply various geometrical transformations of the plane
  to each tile in a tessellation.
}
\usage{
  \method{reflect}{tess}(X)

  \method{shift}{tess}(X, \dots)

  \method{rotate}{tess}(X, angle=pi/2, \dots)

  \method{scalardilate}{tess}(X, f, \dots)

  \method{affine}{tess}(X, mat=diag(c(1,1)), vec=c(0,0), \dots)
}
\arguments{
  \item{X}{Tessellation (object of class \code{"tess"}).}
  \item{angle}{
    Rotation angle in radians (positive values represent
    anticlockwise rotations).
  }
  \item{mat}{Matrix representing a linear transformation.}
  \item{vec}{Vector of length 2 representing a translation.}
  \item{f}{Positive number giving scale factor.}
  \item{\dots}{Arguments passed to other methods.}
}
\value{
  Another tessellation (of class \code{"tess"}) representing the
  result of applying the geometrical transformation.
}
\details{
  These are method for the generic functions 
\code{\link{reflect}},
\code{\link{shift}},
\code{\link{rotate}}, 
\code{\link{scalardilate}},
\code{\link{affine}}
for tessellations (objects of class \code{"tess"}).

The individual tiles of the tessellation, and the window
containing the tessellation, are all subjected to the same
geometrical transformation.

The transformations are performed by the corresponding method
for windows (class \code{"owin"}) or images (class \code{"im"})
depending on the type of tessellation.

If the argument \code{origin} is used in \code{shift.tess}
it is interpreted as applying to the window containing the
tessellation. Then all tiles are shifted by the same vector.
}
\seealso{
  Generic functions
  \code{\link{reflect}},
  \code{\link{shift}},
  \code{\link{rotate}}, 
  \code{\link{scalardilate}},
  \code{\link{affine}}.

  Methods for windows:
  \code{\link{reflect.default}},
  \code{\link{shift.owin}},
  \code{\link{rotate.owin}}, 
  \code{\link{scalardilate.owin}},
  \code{\link{affine.owin}}.

  Methods for images:
  \code{\link{reflect.im}},
  \code{\link{shift.im}},
  \code{\link{rotate.im}}, 
  \code{\link{scalardilate.im}},
  \code{\link{affine.im}}.
}
\examples{
  live <- interactive()
  if(live) {
    H <- hextess(letterR, 0.2)
    plot(H)
    plot(reflect(H))
    plot(rotate(H, pi/3))
  } else H <- hextess(letterR, 0.6)

  # shear transformation
  shear <- matrix(c(1,0,0.6,1),2,2)
  sH <- affine(H, shear)
  if(live) plot(sH)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
