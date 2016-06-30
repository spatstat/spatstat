\name{diameter.linnet}
\alias{circumradius.linnet}
\alias{diameter.linnet}
\title{
  Circumradius and Diameter of a Linear Network
}
\description{
  Compute the circumradius or diameter of a linear network
  measured using the shortest path distance.
}
\usage{
\method{circumradius}{linnet}(x, \dots)

\method{diameter}{linnet}(x)
}
\arguments{
  \item{x}{
    Linear network (object of class \code{"linnet"}).
  }
  \item{\dots}{Ignored.}
}
\details{
  The diameter of a linear network (in the shortest path distance)
  is the maximum value of the shortest-path distance between any
  two points \eqn{u} and \eqn{v} on the network.

  The circumradius of a linear network (in the shortest path distance)
  is the minimum value, over all points \eqn{u} on the network,
  of the maximum shortest-path distance from \eqn{u} to another point
  \eqn{v} on the network.
  
  The functions \code{\link{circumradius}} and
  \code{\link{diameter}} are generic;
  the functions \code{circumradius.linnet} and \code{diameter.linnet}
  are the methods for objects of class \code{linnet}.
}
\value{
  A single numeric value.
}
\author{
  \adrian 
  
}
\seealso{
 \code{\link{circumradius}}, \code{\link{diameter}}, 
 \code{\link{linnet}}
}
\examples{
   diameter(simplenet)
   circumradius(simplenet)
}
\keyword{spatial}
\keyword{math}
