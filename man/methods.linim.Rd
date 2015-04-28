\name{methods.linim}
\Rdversion{1.1}
\alias{methods.linim} %DoNotExport
\alias{as.im.linim}
\alias{print.linim}
\alias{affine.linim}
\alias{scalardilate.linim}
\alias{shift.linim}
\title{
   Methods for Images on a Linear Network
}
\description{
  Methods for the class \code{"linim"} of functions on a linear network.
}
\usage{
  \method{print}{linim}(x, \dots)

  \method{as.im}{linim}(X, \dots)

  \method{shift}{linim}(X, \dots)

  \method{scalardilate}{linim}(X, f, \dots, origin=NULL)

  \method{affine}{linim}(X, mat=diag(c(1,1)), vec=c(0,0), \dots)
}
\arguments{
  \item{X,x}{
    A pixel image on a linear network
    (object of class \code{"linim"}).
  }
  \item{\dots}{
    Extra arguments passed to other methods.
  }
  \item{f}{Numeric. Scalar dilation factor.}
  \item{mat}{Numeric matrix representing the linear transformation.}
  \item{vec}{Numeric vector of length 2 specifying the shift vector.}
  \item{origin}{Character string determining a location
    that will be shifted to the origin. Options are
    \code{"centroid"}, \code{"midpoint"} and \code{"bottomleft"}.
    Partially matched. 
  }
}
\details{
  These are methods for the generic function
  \code{\link{print}}
  and the \pkg{spatstat} generic functions
  \code{\link{as.im}},
  \code{\link{shift}}, 
  \code{\link{scalardilate}} and 
  \code{\link{affine}}.

  An object of class \code{"linfun"} represents a
  pixel image defined on a linear network. 

  The method \code{as.im.linim} extracts the pixel values
  and returns a pixel image of class \code{"im"}.

  The methods \code{shift.linim},
  \code{scalardilate.linim} and \code{affine.linim}
  apply geometric transformations to the pixels and the underlying
  linear network, without changing the pixel values.
}
\value{
  For \code{print.linim} the result is \code{NULL}.

  For \code{as.im.linim} the result is an object of class \code{"im"}.
  
  For the geometric transformations \code{shift.linim},
  \code{scalardilate.linim} and \code{affine.linim}, the result is
  another object of class \code{"linim"}.
}
\examples{
   example(linim)
   X
   shift(X, c(1,1))
   scalardilate(X, 2)
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{math}
