\name{methods.linfun}
\Rdversion{1.1}
\alias{methods.linfun} %DoNotExport
\alias{print.linfun}
\alias{plot.linfun}
\alias{as.linim.linfun}
\alias{as.owin.linfun}
\alias{as.function.linfun}
\title{
   Methods for Functions on Linear Network
}
\description{
  Methods for the class \code{"linfun"} of functions on a linear network.
}
\usage{
  \method{print}{linfun}(x, \dots)

  \method{plot}{linfun}(x, \dots, L=NULL,
                     eps = NULL, dimyx = NULL, xy = NULL,
                     main)

  \method{as.linim}{linfun}(X, L, \dots, eps = NULL, dimyx = NULL, xy = NULL)

  \method{as.owin}{linfun}(W, \dots)

  \method{as.function}{linfun}(x, \dots)
}
\arguments{
  \item{X,x,W}{
    A function on a linear network
    (object of class \code{"linfun"}).
  }
  \item{L}{A linear network}
  \item{eps,dimyx,xy}{
    Arguments passed to \code{\link{as.mask}} to control
    the pixel resolution.
  }
  \item{\dots}{
    Extra arguments passed to \code{\link{plot.im}} or
    \code{\link{print.default}}.
  }
  \item{main}{Main title for plot.}
}
\details{
  These are methods for the generic functions
  \code{\link{plot}}, \code{\link{print}}
  and the \pkg{spatstat} generic functions
  \code{\link{as.owin}},
  \code{\link{as.linim}}.

  An object of class \code{"linfun"} represents a
  mathematical function that could be evaluated at any location
  on a linear network. It is essentially an \R \code{function} with some
  extra attributes.

  The method \code{as.linnet.linfun} extracts the linear network
  on which the function is defined.
  The method \code{as.owin.linfun} extracts the two-dimensional spatial
  window containing the linear network.
  The method \code{as.linim.linfun} converts the function to a pixel
  image on the linear network (an object of class \code{"linim"}).
}
\value{
  For \code{print.linfun} the result is \code{NULL}.

  For \code{plot.linfun} the result is the same as
  for \code{\link{plot.im}}.
  
  For the conversion methods, the result is an object of the
  required type: \code{as.linim.linfun} returns an object of
  class \code{"linim"}, and so on.
}
\examples{
   data(letterR)
   X <- runiflpp(3, simplenet)
   f <- nnfun(X)
   f
   plot(f)
   as.function(f)
   as.owin(f)
   as.linim(f)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
