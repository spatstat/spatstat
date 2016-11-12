\name{ragsMultiHard}
\alias{ragsMultiHard}
\title{
  Alternating Gibbs Sampler for Multitype Hard Core Process
}
\description{
  Generate a realisation of the multitype hard core point process
  using the alternating Gibbs sampler.
}
\usage{
ragsMultiHard(beta, hradii, \dots, types=NULL, bmax = NULL,
              periodic=FALSE, ncycles = 100)
}
\arguments{
  \item{beta}{
    First order trend. A numeric vector, a pixel image,
    a function, a list of functions, or a list of pixel images.
  }
  \item{hradii}{
    Matrix of hard core radii between each pair of types.
    Diagonal entries should be \code{0} or \code{NA}.
  }
  \item{types}{
    Vector of all possible types for the multitype point pattern.
  }
  \item{\dots}{
    Arguments passed to \code{\link{rmpoispp}}
    when generating random points. 
}
  \item{bmax}{
    Optional upper bound on \code{beta}.
  }
  \item{periodic}{
    Logical value indicating whether to measure distances in the
    periodic sense, so that opposite sides of the (rectangular) window
    are treated as identical.
  }
  \item{ncycles}{
    Number of cycles of the sampler to be performed.
  }
}
\details{
  The Alternating Gibbs Sampler for a multitype point process
  is an iterative simulation procedure. Each step of the sampler
  updates the pattern of points of a particular type \code{i},
  by drawing a realisation from the conditional distribution of
  points of type \code{i} given the points of all other types.
  Successive steps of the sampler update the points of type 1, then
  type 2, type 3, and so on. 

  This is an experimental implementation which currently works only
  for multitype hard core processes (see \code{\link{MultiHard}})
  in which there is no interaction between points of the same type,
  and for the area-interaction process (see \code{\link{ragsAreaInter}}).

  The argument \code{beta} gives the first order trend for
  each possible type of point. It may be a single number, a numeric
  vector, a \code{function(x,y)}, a pixel image, a list of functions,
  a \code{function(x,y,m)}, or a list of pixel images. 

  The argument \code{hradii} is the matrix of hard core radii
  between each pair of possible types of points. Two points of types
  \code{i} and \code{j} respectively are forbidden to lie closer than
  a distance \code{hradii[i,j]} apart. The diagonal of this matrix must
  contain \code{NA} or \code{0} values, indicating that there is no hard
  core constraint applying between points of the same type.
}
\value{
  A point pattern (object of class \code{"ppp"}).
}
\author{
  \adrian
}
\seealso{
  \code{\link{rags}},
  \code{\link{ragsAreaInter}}
}
\examples{
  b <- c(30,20)
  h <- 0.05 * matrix(c(0,1,1,0), 2, 2)
  ragsMultiHard(b, h, ncycles=10)
  ragsMultiHard(b, h, ncycles=5, periodic=TRUE)
}
\keyword{spatial}
\keyword{datagen}
