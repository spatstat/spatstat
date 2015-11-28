\name{rpoisppOnLines}
\alias{rpoisppOnLines}
\title{Generate Poisson Point Pattern on Line Segments}
\description{
  Given a line segment pattern, generate a Poisson random point pattern
  on the line segments.
}
\usage{
rpoisppOnLines(lambda, L, lmax = NULL, ..., nsim=1)
}
\arguments{
  \item{lambda}{Intensity of the Poisson process.
    A single number, a \code{function(x,y)}, a pixel image
    (object of class \code{"im"}), or a vector of numbers,
    a list of functions, or a list of images.
  }
  \item{L}{Line segment pattern (object of class \code{"psp"})
    on which the points should be generated.
  }
  \item{lmax}{
    Optional upper bound (for increased computational efficiency).
    A known upper bound for the values of \code{lambda},
    if \code{lambda} is a function or a pixel image.
    That is, \code{lmax} should be a number which is 
    known to be greater than or equal to all values of \code{lambda}. 
  }
  \item{\dots}{Additional arguments passed to \code{lambda} if it is a
    function.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
}
\details{
  This command generates a Poisson point process on the one-dimensional
  system of line segments in \code{L}. The result is a point pattern
  consisting of points lying on the line segments in \code{L}.
  The number of random points falling on any given line segment follows
  a Poisson distribution. The patterns of points on different segments
  are independent.

  The intensity \code{lambda} is the expected number of points
  per unit \bold{length} of line segment. 
  It may be constant, or it may depend on spatial location.

  In order to generate an unmarked Poisson process,
  the argument \code{lambda} may be a single number,
  or a \code{function(x,y)}, or a pixel image (object of class
  \code{"im"}). 

  In order to generate a \emph{marked} Poisson process,
  \code{lambda} may be a numeric vector, a list of functions,
  or a list of images, each entry giving the intensity for
  a different mark value.

  If \code{lambda} is not numeric, then the (Lewis-Shedler)
  rejection method is used. 
  The rejection method requires knowledge of \code{lmax},
  the maximum possible value of \code{lambda}. This should be either
  a single number, or a numeric vector of the same length as
  \code{lambda}. If \code{lmax} is not
  given, it will be computed approximately, by sampling many values of
  \code{lambda}.

  If \code{lmax} is given, then it \bold{must} be larger than
  any possible value of \code{lambda}, otherwise the results of the
  algorithm will be incorrect.
}
\value{
  If \code{nsim = 1}, a point pattern
  (object of class \code{"ppp"}) in the same window
  as \code{L}. If \code{nsim > 1}, a list of such point patterns.
}
\seealso{
  \code{\link{psp}},
  \code{\link{ppp}},
  \code{\link{runifpointOnLines}},
  \code{\link{rpoispp}}
}
\examples{
  live <- interactive()
  L <- psp(runif(10), runif(10), runif(10), runif(10),  window=owin())
  if(live) plot(L, main="")

  # uniform intensity
  Y <- rpoisppOnLines(4, L)
  if(live) plot(Y, add=TRUE, pch="+")

  # uniform MARKED process with types 'a' and 'b'
  Y <- rpoisppOnLines(c(a=4, b=5), L)
  if(live) {
    plot(L, main="")
    plot(Y, add=TRUE, pch="+")
  }

  # intensity is a function
  Y <- rpoisppOnLines(function(x,y){ 10 * x^2}, L, 10)
  if(live) { 
    plot(L, main="")
    plot(Y, add=TRUE, pch="+")
  }

  # intensity is an image
  Z <- as.im(function(x,y){10 * sqrt(x+y)}, unit.square())
  Y <- rpoisppOnLines(Z, L, 15)
  if(live) {
   plot(L, main="")
   plot(Y, add=TRUE, pch="+")
  }
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{datagen}
