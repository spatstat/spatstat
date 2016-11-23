\name{rlpp}
\alias{rlpp}
\title{
  Random Points on a Linear Network
}
\description{
  Generates \eqn{n} independent random points
  on a linear network with a specified probability density.
}
\usage{
  rlpp(n, f, \dots, nsim=1, drop=TRUE)
}
\arguments{
  \item{n}{
    Number of random points to generate.
    A nonnegative integer giving the number of points,
    or an integer vector giving the numbers of points of each type.
  }
  \item{f}{
    Probability density (not necessarily normalised).
    A pixel image on a linear network (object of class \code{"linim"})
    or a function on a linear network (object of class \code{"linfun"}).
    Alternatively, \code{f} can be a list of functions or pixel images,
    giving the densities of points of each type.
  }
  \item{\dots}{
    Additional arguments passed to \code{f} if it is a function
    or a list of functions.
  }
  \item{nsim}{Number of simulated realisations to generate.}
  \item{drop}{
    Logical value indicating what to do when \code{nsim=1}.
    If \code{drop=TRUE} (the default), the result is a point pattern.
    If \code{drop=FALSE}, the result is a list with one entry which is a
    point pattern.
  }
}
\details{
  The linear network \code{L}, on which the points will be generated,
  is determined by the argument \code{f}.
  
  If \code{f} is a function, it is converted to a pixel image
  on the linear network, using any additional function arguments
  \code{\dots}.

  If \code{n} is a single integer and \code{f} is a function or pixel image,
  then independent random points are generated on \code{L} with
  probability density proportional to \code{f}.

  If \code{n} is an integer vector and \code{f} is a list of functions
  or pixel images, where \code{n} and \code{f} have the same length,
  then independent random points of several types are generated on
  \code{L}, with \code{n[i]} points of type \code{i} having probability
  density proportional to \code{f[[i]]}.
}
\value{
  If \code{nsim = 1} and \code{drop=TRUE},
  a point pattern on the linear network,
  i.e.\ an object of class \code{"lpp"}.
  Otherwise, a list of such point patterns.
}
\author{
  \adrian 
}
\seealso{
 \code{\link{runiflpp}}
}
\examples{
  g <- function(x, y, seg, tp) { exp(x + 3*y) }
  f <- linfun(g, simplenet)

  rlpp(20, f)

  plot(rlpp(20, f, nsim=3))
}
\keyword{spatial}
\keyword{datagen}
