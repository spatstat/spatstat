\name{circumradius}
\alias{circumradius}
\alias{circumcentre}
\alias{circumcircle}
\alias{circumradius.owin}
\alias{circumcentre.owin}
\alias{circumcircle.owin}
\alias{circumradius.ppp}
\alias{circumcentre.ppp}
\alias{circumcircle.ppp}
\title{
  Smallest Enclosing Circle
}
\description{
  Find the smallest circle enclosing a spatial window
  or other object. Return its radius, or the location of its centre,
  or the circle itself.
}
\usage{
circumradius(x, \dots)

circumcentre(x, \dots)

circumcircle(x, \dots)

\method{circumradius}{owin}(x, \dots)

\method{circumcentre}{owin}(x, \dots)

\method{circumcircle}{owin}(x, \dots)

\method{circumradius}{ppp}(x, \dots)

\method{circumcentre}{ppp}(x, \dots)

\method{circumcircle}{ppp}(x, \dots)
}
\arguments{
  \item{x}{
    A window (object of class \code{"owin"}), or another spatial object. 
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} to determine the
    pixel resolution for the calculation.
  }
}
\details{
  The circumcircle of a spatial region \eqn{W} is the smallest circle
  that contains \eqn{W}. The circumradius is the radius of this circle,
  and the circumcentre is the centre of the circle.
  
  The functions \code{circumcircle}, \code{circumcentre} and
  \code{circumradius} are generic. There are methods
  for objects of class \code{"owin"}, \code{"ppp"} and \code{"linnet"}.
}
\value{
  The result of \code{circumradius} is a single numeric value.

  The result of \code{circumcentre} is a point pattern containing a
  single point.

  The result of \code{circumcircle} is a window representing the
  circumcircle.
}
\author{
  \adrian 
}
\seealso{
 \code{\link{circumradius.linnet}}
}
\examples{
  circumradius(letterR)

  plot(grow.rectangle(Frame(letterR), 0.2), main="", type="n")
  plot(letterR, add=TRUE, col="grey")
  plot(circumcircle(letterR), add=TRUE, border="green", lwd=2)
  plot(circumcentre(letterR), pch="+", cex=2, col="blue", add=TRUE)

  X <- runifpoint(5)
  plot(X)
  plot(circumcircle(X), add=TRUE)
  plot(circumcentre(X), pch="+", cex=2, col="blue", add=TRUE)
}
\keyword{spatial}
\keyword{math}
