\name{Replace.im}
\alias{[<-.im}
\title{Reset Values in Subset of Image}
\description{
  Reset the values in a subset of a pixel image.
}
\usage{
  \method{[}{im}(x, i, j) <- value
}
\arguments{
  \item{x}{
    A two-dimensional pixel image.
    An object of class \code{"im"}.
  }
  \item{i}{
    Object defining the subregion or subset to be replaced.
    Either a spatial window (an object of class \code{"owin"}), or a
    pixel image with logical values, or a point pattern (an object
    of class \code{"ppp"}), or any type of index that applies to a
    matrix, or something that can be converted to a point pattern
    by \code{\link{as.ppp}} (using the window of \code{x}).
  }
  \item{j}{
    An integer or logical vector serving as the column index
    if matrix indexing is being used.  Ignored if \code{i} is
    appropriate to some sort of replacement \emph{other than}
    matrix indexing.
  }
  \item{value}{
     Vector, matrix, factor or pixel image
     containing the replacement values.
     Short vectors will be recycled.
  }
} 
\value{
  The image \code{x} with the values replaced.
}
\details{
  This function changes some of the pixel values in a
  pixel image. The image \code{x} must be an object of class
  \code{"im"} representing a pixel image defined inside a
  rectangle in two-dimensional space (see \code{\link{im.object}}).

  The subset to be changed is determined by the arguments \code{i,j}
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
  the values of the image inside this window are changed.

  If \code{i} is a point pattern (an object of class
  \code{"ppp"}), then the values of the pixel image at the points of
  this pattern are changed.

  If \code{i} does not satisfy any of the conditions above, then
  the algorithm tries to interpret \code{i,j} as indices for the matrix
  \code{as.matrix(x)}. Either \code{i} or \code{j} may be missing or blank.

  If none of the conditions above are met, and if \code{i} is not
  a matrix, then \code{i} is converted into a point pattern
  by \code{\link{as.ppp}(i, W=Window(x))}.
  Again the values of the pixel image at the points of
  this pattern are changed.
}
\section{Warning}{
  If you have a 2-column matrix containing the \eqn{x,y} coordinates
  of point locations, then to prevent this being interpreted as an
  array index, you should convert it to a \code{data.frame}
  or to a point pattern.
}
\seealso{
  \code{\link{im.object}},
  \code{\link{[.im}},
  \code{\link{[}},
  \code{\link{ppp.object}},
  \code{\link{as.ppp}},
  \code{\link{owin.object}}
}
\examples{
 # make up an image
 X <- setcov(unit.square())
 plot(X)

 # a rectangular subset
 W <- owin(c(0,0.5),c(0.2,0.8))
 X[W] <- 2
 plot(X)

 # a polygonal subset
 data(letterR)
 R <- affine(letterR, diag(c(1,1)/2), c(-2,-0.7))
 X[R] <- 3
 plot(X)

 # a point pattern
 P <- rpoispp(20)
 X[P] <- 10
 plot(X)

 # change pixel value at a specific location
 X[list(x=0.1,y=0.2)] <- 7

 # matrix indexing --- single vector index
 X[1:2570] <- 10
 plot(X)

 # matrix indexing using double indices
 X[1:257,1:10] <- 5
 plot(X)

 # matrix indexing using a matrix of indices
 X[cbind(1:257,1:257)] <- 10
 X[cbind(257:1,1:257)] <- 10
 plot(X)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
