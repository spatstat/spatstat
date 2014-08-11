\name{by.im}
\alias{by.im}
\title{Apply Function to Image Broken Down by Factor}
\description{
  Splits a pixel image into sub-images and applies a function to each
  sub-image.
}
\usage{
\method{by}{im}(data, INDICES, FUN, ...)
}
\arguments{
  \item{data}{A pixel image (object of class \code{"im"}).}
  \item{INDICES}{Grouping variable.
    Either a tessellation (object of class \code{"tess"}) or
    a factor-valued pixel image.
  }
  \item{FUN}{Function to be applied to each sub-image of \code{data}.}
  \item{\dots}{Extra arguments passed to \code{FUN}.}
}
\details{
  This is a method for the generic function \code{\link{by}} for
  pixel images (class \code{"im"}).

  The pixel image \code{data} is first divided into sub-images according
  to \code{INDICES}. Then the function \code{FUN} is applied to each subset.
  The results of each computation are returned in a list.

  The grouping variable \code{INDICES} may be either
  \itemize{
    \item a tessellation (object of class \code{"tess"}). Each tile of
    the tessellation delineates a subset of the spatial domain.
    \item a pixel image (object of class \code{"im"}) with factor
    values. The levels of the factor determine subsets of the spatial
    domain.
  }
}
\value{
  A list containing the results of each evaluation of \code{FUN}.
}
\seealso{
  \code{\link{split.im}},
  \code{\link{tess}},
  \code{\link{im}}
}
\examples{
  W <- square(1)
  X <- as.im(function(x,y){sqrt(x^2+y^2)}, W)
  Y <- dirichlet(runifpoint(12, W))
  # mean pixel value in each subset
  unlist(by(X, Y, mean))
  # trimmed mean
  unlist(by(X, Y, mean, trim=0.05))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{manip}
