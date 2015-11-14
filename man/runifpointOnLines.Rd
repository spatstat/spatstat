\name{runifpointOnLines}
\alias{runifpointOnLines}
\title{Generate N Uniform Random Points On Line Segments}
\description{
  Given a line segment pattern, generate a random point pattern
  consisting of \code{n} points uniformly distributed on the
  line segments.
}
\usage{
runifpointOnLines(n, L, nsim=1)
}
\arguments{
  \item{n}{Number of points to generate.}
  \item{L}{Line segment pattern (object of class \code{"psp"})
    on which the points should lie.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
}
\details{
  This command generates a point pattern consisting of \code{n}
  independent random points, each point uniformly distributed on
  the line segment pattern. This means that, for each random point,
  \itemize{
    \item the
    probability of falling on a particular segment is proportional to the
    length of the segment; and
    \item given that the point falls on a particular segment,
    it has uniform probability density along that segment.
  }
  If \code{n} is a single integer, the result is an unmarked point
  pattern containing \code{n} points.
  If \code{n} is a vector of integers, the result is a marked point
  pattern, with \code{m} different types of points, where
  \code{m = length(n)},
  in which there are \code{n[j]} points of type \code{j}.
}
\value{
  If \code{nsim = 1}, 
  a point pattern (object of class \code{"ppp"}) with the same
  window as \code{L}.
  If \code{nsim > 1}, a list of point patterns.
}
\seealso{
  \code{\link{psp}},
  \code{\link{ppp}},
  \code{\link{pointsOnLines}},
  \code{\link{runifpoint}}
}
\examples{
  X <- psp(runif(10), runif(10), runif(10), runif(10),  window=owin())
  Y <- runifpointOnLines(20, X)
  plot(X, main="")
  plot(Y, add=TRUE)
  Z <- runifpointOnLines(c(5,5), X)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{datagen}
