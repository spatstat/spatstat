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

  \code{f} should be a \code{function} in the \R language,
  with formal arguments \code{f(x,y,seg,tp)} or 
  \code{f(x,y,seg,tp, \dots)} where \code{x,y} are
  Cartesian coordinates of locations on the linear network,
  \code{seg, tp} are the local coordinates, and
  \code{\dots} are optional additional arguments.
  
  The function \code{f} should be vectorised: that is,
  if \code{x,y,seg,tp} are numeric vectors of the same length
  \code{n}, then \code{v <- f(x,y,seg,tp)}
  should be a vector of length \code{n}.

  \code{L} should be a linear network (object of class \code{"linnet"})
  inside which the function \code{f} is well-defined.

  The result is a function \code{g} in the \R language which belongs to
  the special class \code{"linfun"}. This function
  can be called as \code{g(X)} where \code{X} is an \code{"lpp"} object,
  or called as \code{g(x,y)} or \code{g(x,y,seg,tp)} where
  \code{x,y,seg,tp} are coordinates. There are several methods
  for this class including \code{print}, \code{plot}
  and \code{\link{as.linim}}.
}
\value{
  A function in the \R\ language.
  It also belongs to the class \code{"linfun"} which has methods
  for \code{plot}, \code{print} etc.
}
\seealso{
  \code{\link{methods.linfun}} for methods applicable to
  \code{"linfun"} objects.

  \code{\link{distfun.lpp}},
  \code{\link{nnfun.lpp}}.
}
\examples{
  f <- linfun(function(x,y,seg,tp) { x+y }, simplenet)
  plot(f)
  X <- runiflpp(3, simplenet)
  plot(X, add=TRUE, cex=2)
  f(X)
}
\author{\adrian
  and \rolf
}
\keyword{spatial}
\keyword{math}
