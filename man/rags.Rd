\name{rags}
\alias{rags}
\title{
  Alternating Gibbs Sampler for Multitype Point Processes
}
\description{
  Simulate a realisation of a point process model using the
  alternating Gibbs sampler.
}
\usage{
rags(model, \dots, ncycles = 100)
}
\arguments{
  \item{model}{
    Data specifying some kind of point process model.
  }
  \item{\dots}{
    Additional arguments passed to other code.
  }
  \item{ncycles}{
    Number of cycles of the alternating Gibbs sampler that should be
    performed. 
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
  in which there is no interaction between points of the same type. 

  The argument \code{model} should be an object describing a point
  process model. At the moment, the only permitted format for
  \code{model} is of the form \code{list(beta, hradii)} where
  \code{beta} gives the first order trend and \code{hradii} is the
  matrix of interaction radii. See \code{\link{ragsMultiHard}} for
  full details.
}
\value{
  A point pattern (object of class \code{"ppp"}).
}
\author{
  \adrian
}
\seealso{
  \code{\link{ragsMultiHard}},
  \code{\link{ragsAreaInter}}
}
\examples{
  mo <- list(beta=c(30, 20),
             hradii = 0.05 * matrix(c(0,1,1,0), 2, 2))
  rags(mo, ncycles=10)
}
\keyword{spatial}
\keyword{datagen}
