\name{nnfun}
\Rdversion{1.1}
\alias{nnfun}
\alias{nnfun.ppp}
\alias{nnfun.psp}
\title{
  Nearest Neighbour Index Map as a Function 
}
\description{
  Compute the nearest neighbour index map of an object,
  and return it as a function.
}
\usage{
  nnfun(X, ...)

  \method{nnfun}{ppp}(X, ..., k=1)

  \method{nnfun}{psp}(X, ...)
}
\arguments{
  \item{X}{
    Any suitable dataset representing a two-dimensional
    collection of objects,
    such as a point pattern (object of class \code{"ppp"})
    or a line segment pattern (object of class \code{"psp"}).
  }
  \item{k}{
    A single integer. The \code{k}th nearest neighbour will be found.
  }
  \item{\dots}{
    Extra arguments are ignored. 
  }
}
\details{
  For a collection \eqn{X} of two dimensional objects
  (such as a point pattern or a line segment pattern),
  the \dQuote{nearest neighbour index function}
  of \eqn{X} is the mathematical function \eqn{f} such that, for any 
  two-dimensional spatial location \eqn{(x,y)},
  the function value \code{f(x,y)}
  is the index \eqn{i} identifying the closest member of \eqn{X}.
  That is, if \eqn{i = f(x,y)} then \eqn{X[i]} is the closest member of
  the collection \eqn{X} to the location \eqn{(x,y)}.
  
  The command \code{f <- nnfun(X)} returns a \emph{function}
  in the \R language, with arguments \code{x,y}, that represents the
  nearest neighbour index function of \code{X}. Evaluating the function \code{f}
  in the form \code{v <- f(x,y)}, where \code{x} and \code{y}
  are any numeric vectors of equal length containing coordinates of
  spatial locations, yields the indices of the nearest neighbours
  to these locations.

  If the argument \code{k} is specified then the \code{k}-th nearest
  neighbour will be found.

  The result of \code{f <- nnfun(X)} also belongs to the class
  \code{"funxy"} and to the special class \code{"nnfun"}.
  It can be printed and plotted immediately as shown in the Examples.
  
  A \code{nnfun} object can be converted to a pixel image
  using \code{\link{as.im}}.
}
\value{
  A \code{function} with arguments \code{x,y}.
  The function also belongs to the class \code{"nnfun"} which has
  a method for \code{print}.
  It also belongs to the class \code{"funxy"} which has methods
  for \code{plot}, \code{contour} and \code{persp}.
}
\seealso{
  \code{\link{distfun}},
  \code{\link{plot.funxy}}
}
\examples{
   f <- nnfun(cells)
   f
   plot(f)
   f(0.2, 0.3)

   g <- nnfun(cells, k=2)
   g(0.2, 0.3)

   L <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
   h <- nnfun(L)
   h(0.2, 0.3)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{math}
