\name{polartess}
\alias{polartess}
\title{
  Tessellation Using Polar Coordinates
}
\description{
  Create a tessellation with tiles defined by polar coordinates
  (radius and angle).
}
\usage{
polartess(W, \dots, nradial = NULL, nangular = NULL,
          radii = NULL, angles = NULL,
          origin = NULL, sep = "x")
}
\arguments{
  \item{W}{
    A window (object of class \code{"owin"})
    or anything that can be coerced to a window using
    \code{\link{as.owin}}, such as a point pattern.
  }
  \item{\dots}{
    Ignored.
  }
  \item{nradial}{
    Number of \emph{tiles} in the radial direction.
    A single integer.
    Ignored if \code{radii} is given.
  }
  \item{nangular}{
    Number of \emph{tiles} in the angular coordinate.
    A single integer.
    Ignored if \code{angles} is given.
  }
  \item{radii}{
    The numeric values of the radii, defining the tiles in the radial
    direction. A numeric vector, of length at least 2, containing
    nonnegative numbers in increasing order. The value \code{Inf} is permitted.
  }
  \item{angles}{
    The numeric values of the angles defining the tiles in the
    angular coordinate. A numeric vector, of length at least 2,
    in increasing order, containing angles in radians.
  }
  \item{origin}{
    Location to be used as the origin of the polar coordinates.
    Either a numeric vector of length 2 giving the spatial location
    of the origin, or one of the strings
    \code{"centroid"}, \code{"midpoint"},
    \code{"left"}, \code{"right"}, \code{"top"}, \code{"bottom"},
    \code{"topleft"}, \code{"bottomleft"}, \code{"topright"} or
    \code{"bottomright"} indicating the location in the window.
  }
  \item{sep}{
    Argument passed to \code{\link{intersect.tess}}
    specifying the character string to be used as a separator
    when forming the names of the tiles.
  }
}
\details{
  A tessellation will be formed from tiles defined by
  intervals in the polar coordinates \eqn{r} (radial distance
  from the origin) or \eqn{\theta} (angle from the horizontal axis)
  or both. These tiles look like the cells on a dartboard.

  If the argument \code{radii} is given, tiles will be demarcated
  by circles centred at the origin, with the specified radii.
  If \code{radii} is absent but \code{nradial} is given,
  then \code{radii} will default to a sequence
  of \code{nradial+1} radii equally spaced from zero to the maximum
  possible radius. If neither \code{radii} nor \code{nradial} are given,
  the tessellation will not include circular arc boundaries.

  If the argument \code{angles} is given, tiles will be demarcated
  by lines emanating from the origin at the specified angles.
  The angular values can be any real numbers; they will be interpreted
  as angles in radians modulo \code{2*pi},
  but they must be an increasing sequence of numbers.
  If \code{angles} is absent but \code{nangular} is given,
  then \code{angles} will default to a sequence
  of \code{nangular+1} angles equally spaced from 0 to \code{2*pi}.
  If neither \code{angles} nor \code{nangular} are given,
  the tessellation will not include linear boundaries.
}
\value{
  A tessellation (object of class \code{"tess"}).
}
\author{
  \adrian.
}
\seealso{
  \code{\link{intersect.tess}}
  
  To construct other kinds of tessellations, see
  \code{\link{tess}},
  \code{\link{quadrats}}, \code{\link{hextess}},
  \code{\link{venn.tess}},
  \code{\link{dirichlet}}, \code{\link{delaunay}}, \code{\link{quantess}}
  and \code{\link{rpoislinetess}}.
  
}
\examples{
  Y <- c(2.8, 1.5)
  plot(polartess(letterR, nangular=6, radii=(0:4)/2, origin=Y),
       do.col=TRUE)
}
\keyword{spatial}
\keyword{manip}
