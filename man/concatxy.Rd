\name{concatxy}
\alias{concatxy}
\title{Concatenate x,y Coordinate Vectors}
\description{
  Concatenate any number of pairs of \code{x} and \code{y}
  coordinate vectors.
}
\usage{
  concatxy(\dots)
}
\arguments{
  \item{\dots}{
    Any number of arguments, each of which is a structure
    containing elements \code{x} and \code{y}.
  }
}
\value{
  A list with two components \code{x} and \code{y}, which are the
  concatenations of all the corresponding \code{x} and \code{y}
  vectors in the argument list.
}
\details{
  This function can be used to superimpose two or more point patterns
  of unmarked points (but see also \code{\link{superimpose}} which is
  recommended). 
  
  It assumes that each of the arguments in
  \code{\dots} is a structure containing (at least) the elements
  \code{x} and \code{y}. It concatenates all the \code{x} elements
  into a vector \code{x}, and similarly for \code{y}, and returns these
  concatenated vectors.
}
\seealso{
  \code{\link{superimpose}},
  \code{\link{quadscheme}}
}
\examples{
  dat <- runifrect(30)
  xy <- list(x=runif(10),y=runif(10))
  new <- concatxy(dat, xy)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
