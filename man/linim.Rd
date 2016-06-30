\name{linim}
\alias{linim}
\title{
  Create Pixel Image on Linear Network
}
\description{
  Creates an object of class \code{"linim"} that represents
  a pixel image on a linear network.
}
\usage{
  linim(L, Z, ..., df=NULL)
}
\arguments{
  \item{L}{
    Linear network (object of class \code{"linnet"}).
  }
  \item{Z}{
    Pixel image (object of class \code{"im"}).
  }
  \item{\dots}{Ignored.}
  \item{df}{
    Advanced use only. Data frame giving full details of the mapping between
    the pixels of \code{Z} and the lines of \code{L}.
    See Details.
  }
}
\details{
  This command creates an object of class \code{"linim"} that represents
  a pixel image defined on a linear network.
  Typically such objects are
  used to represent the result of smoothing or model-fitting on the
  network. Most users will not need to call \code{linim} directly.

  The argument \code{L} is a linear network (object of class \code{"linnet"}).
  It gives the exact spatial locations
  of the line segments of the network, and their connectivity.

  The argument \code{Z} is a pixel image object of class \code{"im"}
  that gives a pixellated approximation of the function values.
  
  For increased efficiency, advanced users may specify the 
  optional argument \code{df}. This is a data frame giving the
  precomputed mapping between the pixels of \code{Z}
  and the line segments of \code{L}.
  It should have columns named \code{xc, yc} containing the coordinates of
  the pixel centres, \code{x,y} containing the projections of these
  pixel centres onto the linear network, \code{mapXY} identifying the
  line segment on which each projected point lies, and \code{tp} giving
  the parametric position of \code{(x,y)} along the segment.
}
\value{
  Object of class \code{"linim"} that also inherits the class
  \code{"im"}.
  There is a special method for plotting this class.
}
\author{
  \adrian
  
  
}
\seealso{
  \code{\link{plot.linim}},
  \code{\link{linnet}},
  \code{\link{eval.linim}},
  \code{\link{Math.linim}},
  \code{\link{im}}.
}
\examples{
  M <- as.mask.psp(as.psp(simplenet))
  Z <- as.im(function(x,y) {x-y}, W=M)
  X <- linim(simplenet, Z)
  X
}
\references{
  Ang, Q.W. (2010)
  \emph{Statistical methodology for events on a network}.
  Master's thesis, School of Mathematics and Statistics, University of
  Western Australia.
  
  Ang, Q.W., Baddeley, A. and Nair, G. (2012)
  Geometrically corrected second-order analysis of 
  events on a linear network, with applications to
  ecology and criminology.
  \emph{Scandinavian Journal of Statistics} \bold{39}, 591--617.

  McSwiggan, G., Nair, M.G. and Baddeley, A. (2012)
  Fitting Poisson point process models to events 
  on a linear network. Manuscript in preparation.
}
\keyword{spatial}
