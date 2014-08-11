\name{linfun}
\Rdversion{1.1}
\alias{linfun}
\title{
  Function on a Linear Network
}
\description{
  Create a function on a linear network.
}
\usage{
  linfun(f, L)
}
\arguments{
  \item{f}{
    A \code{function} in the \R language.
  }
  \item{L}{
    A linear network (object of class \code{"linnet"})
    on which \code{f} is defined.
  }
}
\details{
  This creates an object of class \code{"linfun"}.
  This is a simple mechanism for handling a function
  defined on a linear network, to make it easier to display
  and manipulate.

  \code{f} should be a \code{function} in the \R language.
  It should be able to be called either in the form
  \code{f(x,y)} or \code{f(x,y,seg,tp)} where \code{x,y} are
  Cartesian coordinates of locations on the linear network,
  and \code{seg, tp} are the local coordinates.
  The function \code{f} should be vectorised: that is,
  if \code{x} and \code{y} are numeric vectors of the same length
  \code{n}, then \code{v <- f(x,y)} should be a vector of length \code{n}.

  \code{L} should be a linear network (object of class \code{"linnet"})
  inside which the function \code{f} is well-defined.
}
\value{
  It also belongs to the class \code{"linfun"} which has methods
  for \code{plot}, \code{print} etc.
}
\seealso{
  \code{\link{as.lpp}} for a way to handle the coordinates
  \code{x,y,seg,tp}.
  
  \code{\link{methods.linfun}} for methods applicable to
  \code{"linfun"} objects.

  \code{\link{distfun.lpp}},
  \code{\link{nnfun.lpp}}.
}
\examples{
  f <- linfun(function(x,y,seg,tp) { x+y }, simplenet)
  plot(f)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
