\name{discs}
\alias{discs}
\title{
 Union of Discs
}
\description{
  Make a spatial region composed of discs
  with given centres and radii.  
}
\usage{
  discs(centres, radii = marks(centres)/2, \dots,
        separate = FALSE, mask = FALSE, trim = TRUE, delta = NULL)
}
\arguments{
  \item{centres}{
    Point pattern giving the locations of centres for the discs.
  }
  \item{radii}{
    Vector of radii for each disc, or a single number giving a common
    radius.
    (Notice that the default assumes that the marks of \code{X} are
    \emph{diameters}.)
  }
  \item{\dots}{
    Optional arguments passed to \code{\link{as.mask}} to determine
    the pixel resolution, if \code{mask=TRUE}.
  }
  \item{separate}{
    Logical. If \code{TRUE}, the result is a list
    containing each disc as a separate entry.
    If \code{FALSE} (the default), the result is a window
    obtained by forming the union of the discs.
  }
  \item{mask}{
    Logical. If \code{TRUE}, the result is a binary mask window.
    If \code{FALSE}, the result is a polygonal window.
    Applies only when \code{separate=FALSE}.
  }
  \item{trim}{
    Logical value indicating whether to restrict the result
    to the original window of the \code{centres}.
    Applies only when \code{separate=FALSE}.
  }
  \item{delta}{
    Argument passed to \code{\link{disc}} to determine the 
    tolerance for the polygonal approximation of each disc.
    Applies only when \code{mask=FALSE}.
  }
}
\details{
  This command is typically applied to a marked point pattern
  dataset \code{X} in which the marks represent the sizes of objects.
  The result is a spatial region representing the space occupied by
  the objects.
  
  If the marks of \code{X} represent the diameters of circular objects,
  then the result of \code{discs(X)}
  is a spatial region constructed by taking discs, of the specified
  diameters, centred at the points of \code{X}, and forming the union
  of these discs. If the marks of \code{X} represent the areas of
  objects, one could take \code{discs(X, sqrt(marks(X)/pi))}
  to produce discs of equivalent area.

  A fast algorithm is used to compute the result as a binary mask, when
  \code{mask=TRUE}. This option is recommended unless polygons are
  really necessary.

  If \code{mask=FALSE}, the discs will be constructed as polygons
  by the function \code{\link{disc}}. To avoid computational problems, the discs
  will all be constructed using the same physical tolerance value \code{delta}
  passed to \code{\link{disc}}. The default is such that the smallest
  disc will be approximated by a 16-sided polygon.
}
\value{
  If \code{separate=FALSE}, a window (object of class \code{"owin"}).

  If \code{separate=TRUE}, a list of windows.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\seealso{
  \code{\link{disc}},
  \code{\link{union.owin}}
}
\examples{
  plot(discs(anemones, mask=TRUE, eps=0.5))
}
\keyword{spatial}
\keyword{datagen}
