\name{rstrat}
\alias{rstrat}
\title{Simulate Stratified Random Point Pattern}
\description{
  Generates a ``stratified random'' pattern of points in a window,
  by dividing the window into rectangular tiles and placing
  \code{k} random points independently in each tile.
}
\usage{
 rstrat(win=square(1), nx, ny=nx, k = 1, nsim=1, drop=TRUE)
}
\arguments{
  \item{win}{
    A window. 
    An object of class \code{\link{owin}},
    or data in any format acceptable to \code{\link{as.owin}()}.
  }
  \item{nx}{Number of tiles in each column.
  }
  \item{ny}{Number of tiles in each row.
  }
  \item{k}{Number of random points to generate in each tile.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\value{
  A point pattern (an object of class \code{"ppp"})
  if \code{nsim=1}, or a list of point patterns if \code{nsim > 1}.
}
\details{
  This function generates a random pattern of points
  in a ``stratified random'' sampling design. It can be useful
  for generating random spatial sampling points.

  The bounding rectangle of \code{win} is divided into
  a regular \eqn{nx \times ny}{nx * ny} grid of rectangular tiles.
  In each tile, \code{k} random points are generated independently
  with a uniform distribution in that tile. 

  Some of these grid points may lie outside the window \code{win}:
  if they do, they are deleted.

  The result is a point pattern inside the window \code{win}.

  This function is useful in creating dummy points for quadrature
  schemes (see \code{\link{quadscheme}}) as well as in simulating
  random point patterns.
}
\seealso{
  \code{\link{rsyst}},
  \code{\link{runifpoint}},
  \code{\link{quadscheme}}
}
\examples{
  X <- rstrat(nx=10)
  plot(X)

  # polygonal boundary
  data(letterR)
  X <- rstrat(letterR, 5, 10, k=3)
  plot(X)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}

