\name{nnfun.lpp}
\Rdversion{1.1}
\alias{nnfun.lpp}
\title{
   Nearest Neighbour Map on Linear Network
}
\description{
  Compute the nearest neighbour function of a point pattern on a linear network.
}
\usage{
  \method{nnfun}{lpp}(X, ...)
}
\arguments{
  \item{X}{
    A point pattern on a linear network
    (object of class \code{"lpp"}).
  }
  \item{\dots}{
    Extra arguments are ignored. 
  }
}
\details{
  The (geodesic) \emph{nearest neighbour function} of a
  point pattern \code{X} on a linear network \code{L}
  tells us which point of \code{X} is closest to
  any given location.
  
  If \code{X} is a point pattern on a linear network \code{L},
  the \emph{nearest neighbour function} of \code{X}
  is the mathematical function \eqn{f} defined for any 
  location \eqn{s} on the network by \code{f(s) = i}, where
  \code{X[i]} is the closest point of \code{X} to the location \code{s}
  measured by the shortest path. In other words the value of \code{f(s)}
  is the identifier or serial number of the closest point of \code{X}.
  
  The command \code{nnfun.lpp} is a method for the generic command
  \code{\link{nnfun}}
  for the class \code{"lpp"} of point patterns on a linear network.

  If \code{X} is a point pattern on a linear network,
  \code{f <- nnfun(X)} returns a \emph{function}
  in the \R language, with arguments \code{x,y, \dots}, that represents the
  nearest neighbour function of \code{X}. Evaluating the function \code{f}
  in the form \code{v <- f(x,y)}, where \code{x} and \code{y}
  are any numeric vectors of equal length containing coordinates of
  spatial locations, yields a vector of identifiers or serial numbers of
  the data points closest to these spatial locations.
  More efficiently \code{f} can take the arguments
  \code{x, y, seg, tp} where \code{seg} and \code{tp} are the local
  coordinates on the network.

  The result of \code{f <- nnfun(X)} also belongs to the class
  \code{"linfun"}.
  It can be printed and plotted immediately as shown in the Examples.
  It can be converted to a pixel image
  using \code{\link{as.linim}}.

}
\value{
  A \code{function} in the \R language, with arguments \code{x,y} and optional
  arguments \code{seg,tp}.
  It also belongs to the class \code{"linfun"} which has methods
  for \code{plot}, \code{print} etc.
}
\seealso{
  \code{\link{linfun}},
  \code{\link{methods.linfun}}.

  To compute the \emph{distance} to the nearest neighbour, see
  \code{\link{distfun.lpp}}.
}
\examples{
   data(letterR)
   X <- runiflpp(3, simplenet)
   f <- nnfun(X)
   f
   plot(f)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
