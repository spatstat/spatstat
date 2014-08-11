\name{circumradius}
\alias{circumradius}
\alias{diameter.linnet}
\title{
  Circumradius and Diameter of a Linear Network
}
\description{
  Compute the circumradius or diameter of a linear network
  measured using the shortest path distance.
}
\usage{
circumradius(x)
\method{diameter}{linnet}(x)
}
\arguments{
  \item{x}{
    Linear network (object of class \code{"linnet"}).
  }
}
\details{
  The diameter of a linear network (in the shortest path distance)
  is the maximum value of the shortest-path distance between any
  two points \eqn{u} and \eqn{v} on the network.

  The circumradius of a linear network (in the shortest path distance)
  is the minimum value, over all points \eqn{u} on the network,
  of the maximum shortest-path distance from \eqn{u} to another point
  \eqn{v} on the network.
  
  The function \code{diameter} is generic; the function
  \code{diameter.linnet} is the method for objects of class \code{linnet}.
}
\value{
  A single numeric value.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
 \code{\link{linnet}}
}
\examples{
   data(simplenet)
   diameter(simplenet)
   circumradius(simplenet)
}
\keyword{spatial}
\keyword{math}
