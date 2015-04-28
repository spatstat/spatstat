\name{dirichletWeights}
\alias{dirichletWeights}
\title{Compute Quadrature Weights Based on Dirichlet Tessellation}
\description{
  Computes quadrature weights for a given set of points,
  using the areas of tiles in the Dirichlet tessellation.
}
\usage{
 dirichletWeights(X, window=NULL, exact=TRUE, \dots)
}
\arguments{
  \item{X}{Data defining a point pattern.}
  \item{window}{Default window for the point pattern}
  \item{exact}{Logical value. If \code{TRUE}, compute exact areas
    using the package \code{deldir}. If \code{FALSE}, compute
    approximate areas using a pixel raster.
  }
  \item{\dots}{
    Ignored.
  }
}
\value{
  Vector of nonnegative weights for each point in \code{X}.
}
\details{
  This function computes a set of quadrature weights
  for a given pattern of points
  (typically comprising both ``data'' and `dummy'' points).
  See \code{\link{quad.object}} for an explanation of quadrature
  weights and quadrature schemes.

  The weights are computed using the Dirichlet tessellation.
  First \code{X} and (optionally) \code{window} are converted into a
  point pattern object. Then the Dirichlet tessellation of the points
  of \code{X} is computed.
  The weight attached to a point of \code{X} is the area of
  its Dirichlet tile (inside the window \code{X$window}).

  If \code{exact=TRUE} the Dirichlet tessellation is computed exactly
  by the Lee-Schachter algorithm using the package \code{deldir}.
  Otherwise a pixel raster approximation is constructed and the areas
  are approximations to the true weights. In all cases the sum of the
  weights is equal to the area of the window.
}
\seealso{
  \code{\link{quad.object}},
  \code{\link{gridweights}}
}
\examples{
  Q <- quadscheme(runifpoispp(10))
  X <- as.ppp(Q) # data and dummy points together
  w <- dirichletWeights(X, exact=FALSE)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{utilities}
