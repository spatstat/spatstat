\name{split.im}
\alias{split.im}
\title{Divide Image Into Sub-images}
\description{
  Divides a pixel image into several sub-images according to the
  value of a factor, or according to the tiles of a tessellation.
}
\usage{
\method{split}{im}(x, f, ..., drop = FALSE)
}
\arguments{
  \item{x}{Pixel image (object of class \code{"im"}).}
  \item{f}{
    Splitting criterion. Either a tessellation (object of class
    \code{"tess"}) or a pixel image with factor values.
  }
  \item{\dots}{Ignored.}
  \item{drop}{Logical value determining whether each subset should
    be returned as a pixel images (\code{drop=FALSE}) or 
    as a one-dimensional vector of pixel values (\code{drop=TRUE}).
  }
}
\details{
  This is a method for the generic function \code{\link{split}}
  for the class of pixel images. The image \code{x} will be divided
  into subsets determined by the data \code{f}. The result is a list
  of these subsets.

  The splitting criterion may be either
  \itemize{
    \item a tessellation (object of class \code{"tess"}). Each tile of
    the tessellation delineates a subset of the spatial domain.
    \item a pixel image (object of class \code{"im"}) with factor
    values. The levels of the factor determine subsets of the spatial
    domain.
  }
  If \code{drop=FALSE} (the default), the result is a list of pixel
  images, each one a subset of the pixel image \code{x},
  obtained by restricting the pixel domain to one of the subsets.
  If \code{drop=TRUE}, then the pixel values are returned as
  numeric vectors. 
}
\value{
  If \code{drop=FALSE}, a list of pixel images (objects of class
  \code{"im"}). It is also of class \code{"listof"} so that it can be
  plotted immediately.

  If \code{drop=TRUE}, a list of numeric vectors.
}
\seealso{
  \code{\link{by.im}},
  \code{\link{tess}},
  \code{\link{im}}
}
\examples{
  W <- square(1)
  X <- as.im(function(x,y){sqrt(x^2+y^2)}, W)
  Y <- dirichlet(runifpoint(12, W))
  plot(split(X,Y))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{manip}
