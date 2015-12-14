\name{convexify}
\alias{convexify}
\title{
  Weil's Convexifying Operation
}
\description{
  Converts the window \code{W} into a convex set by rearranging
  the edges, preserving spatial orientation of each edge.
}
\usage{
 convexify(W, eps)
}
\arguments{
  \item{W}{
    A window (object of class \code{"owin"}).
  }
  \item{eps}{
    Optional. Minimum edge length of polygonal approximation,
    if \code{W} is not a polygon.
  }
}
\details{
  Weil (1995) defined a convexification operation 
  for windows \eqn{W} that belong to the convex ring (that is,
  for any \eqn{W} which is a finite union of convex sets).
  Note that this is \bold{not} the same as the convex hull.

  The convexified set \eqn{f(W)} has the same total boundary length as
  \eqn{W} and the same distribution of orientations of the boundary.
  If \eqn{W} is a polygonal set, then the convexification \eqn{f(W)}
  is obtained by rearranging all the edges of \eqn{W} in order of
  their spatial orientation.

  The argument \code{W} must be a window. If it is not already a polygonal
  window, it is first converted to one, using
  \code{\link{simplify.owin}}.
  The edges are sorted in increasing order of angular orientation
  and reassembled into a convex polygon.
}
\value{
  A window (object of class \code{"owin"}).
}
\references{
  Weil, W. (1995)
  The estimation of mean particle shape and mean
  particle number in overlapping particle systems in
  the plane. \emph{Advances in Applied Probability} \bold{27}, 102--119.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\seealso{
  \code{\link{convexhull}} for the convex hull of a window.
}
\examples{
  opa <- par(mfrow=c(1,2))
  plot(letterR)
  plot(convexify(letterR))
  par(opa)
}
\keyword{spatial}
\keyword{utilities}
