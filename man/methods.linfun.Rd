\name{methods.linfun}
\Rdversion{1.1}
\alias{methods.linfun} %DoNotExport
\alias{print.linfun}
\alias{summary.linfun}
\alias{plot.linfun}
\alias{as.data.frame.linfun}
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

  \method{summary}{linfun}(object, \dots)

  \method{plot}{linfun}(x, \dots, L=NULL, main) 

  \method{as.data.frame}{linfun}(x, \dots)

  \method{as.owin}{linfun}(W, \dots)

  \method{as.function}{linfun}(x, \dots)
}
\arguments{
  \item{x,object,W}{
    A function on a linear network
    (object of class \code{"linfun"}).
  }
  \item{L}{A linear network}
  \item{\dots}{
    Extra arguments passed to \code{\link{as.linim}}, 
    \code{\link{plot.linim}}, \code{\link{plot.im}}
    or \code{\link{print.default}}, or arguments passed to
    \code{x} if it is a function.
  }
  \item{main}{Main title for plot.}
}
\details{
  These are methods for the generic functions
  \code{\link{plot}}, \code{\link{print}}, \code{\link{summary}}
  \code{\link{as.data.frame}} and \code{\link{as.function}},
  and for the \pkg{spatstat} generic function
  \code{\link{as.owin}}.

  An object of class \code{"linfun"} represents a
  mathematical function that could be evaluated at any location
  on a linear network. It is essentially an \R \code{function} with some
  extra attributes.

  The method \code{as.owin.linfun} extracts the two-dimensional spatial
  window containing the linear network.

  The method \code{plot.linfun} first converts the function to a
  pixel image using \code{\link{as.linim.linfun}}, then plots the image using
  \code{\link{plot.linim}}.

  Note that a \code{linfun} function may have additional arguments,
  other than those which specify the location on the network
  (see \code{\link{linfun}}). These additional arguments may be passed
  to \code{plot.linfun}. 
}
\value{
  For \code{print.linfun} and \code{summary.linfun} the result is \code{NULL}.

  For \code{plot.linfun} the result is the same as
  for \code{\link{plot.im}}.
  
  For the conversion methods, the result is an object of the
  required type: \code{as.owin.linfun} returns an object of
  class \code{"owin"}, and so on.
}
\examples{
   X <- runiflpp(3, simplenet)
   f <- nnfun(X)
   f
   plot(f)
   as.function(f)
   as.owin(f)
   head(as.data.frame(f))
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{math}
