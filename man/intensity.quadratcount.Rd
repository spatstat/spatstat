\name{intensity.quadratcount} 
\alias{intensity.quadratcount}
\title{
  Intensity Estimates Using Quadrat Counts
}
\description{
  Uses quadrat count data to estimate the intensity of a point pattern
  in each tile of a tessellation, assuming the intensity is constant in
  each tile.
}
\usage{
\method{intensity}{quadratcount}(X, ..., image=FALSE)
}
\arguments{
  \item{X}{
    An object of class \code{"quadratcount"}.
  }
  \item{image}{
    Logical value specifying whether to return
    a table of estimated intensities (the default)
    or a pixel image of the estimated intensity (\code{image=TRUE}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} to determine
    the resolution of the pixel image, if \code{image=TRUE}.
  }
}
\details{
  This is a method for the generic function \code{\link{intensity}}.
  It computes an estimate of the intensity of a point pattern
  from its quadrat counts.

  The argument \code{X} should be an object of class
  \code{"quadratcount"}. It would have been obtained by applying the function
  \code{\link{quadratcount}} to a point pattern
  (object of class \code{"ppp"}). It contains
  the counts of the numbers of points of the point pattern falling in each
  tile of a tessellation.

  Using this information, \code{intensity.quadratcount}
  divides the quadrat counts by the tile areas,
  yielding the average density of points per unit area
  in each tile of the tessellation.

  If \code{image=FALSE} (the default), these intensity values
  are returned in a contingency table. Cells of the contingency
  table correspond to tiles of the tessellation.

  If \code{image=TRUE}, the estimated intensity function is
  returned as a pixel image. For each pixel, the pixel value is the
  estimated intensity in the tile which contains that pixel.
}
\value{
  If \code{image=FALSE} (the default), a contingency table.
  If \code{image=TRUE}, a pixel image (object of class \code{"im"}).
}
\seealso{
  \code{\link{intensity}}, 
  \code{\link{quadratcount}}
}
\examples{
  qa <- quadratcount(swedishpines, 4,3)
  qa
  intensity(qa)
  plot(intensity(qa, image=TRUE))
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{nonparametric}
