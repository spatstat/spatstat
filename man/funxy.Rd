\name{funxy}
\Rdversion{1.1}
\alias{funxy}
\title{
  Spatial Function Class
}
\description{
  A simple class of functions of spatial location
}
\usage{
  funxy(f, W)
}
\arguments{
  \item{f}{
    A \code{function} in the \R language
    with arguments \code{x,y} (at least)
  }
  \item{W}{
    Window (object of class \code{"owin"}) inside which the
    function is well-defined.
  }
}
\details{
  This creates an object of class \code{"funxy"}.
  This is a simple mechanism for handling a function
  of spatial location \eqn{f(x,y)} to make it easier to display
  and manipulate.

  \code{f} should be a \code{function} in the \R language.
  The first two arguments of \code{f} must be named \code{x} and \code{y}
  respectively.

  \code{W} should be a window (object of class \code{"owin"}) inside which the
  function \code{f} is well-defined.

  The function \code{f} should be vectorised: that is,
  if \code{x} and \code{y} are numeric vectors of the same length
  \code{n}, then \code{v <- f(x,y)} should be a vector of length
  \code{n}.

  The resulting function \code{g <- funxy(f, W)} has the same formal
  arguments as \code{f} and can be called in the same way,
  \code{v <- g(x,y)} where \code{x} and \code{y} are numeric vectors.
  However it can also be called as \code{v <- g(X)}, where \code{X} is a point pattern
  (object of class \code{"ppp"} or \code{"lpp"})
  or a quadrature scheme (class \code{"quad"});
  the function will be evaluated at the points of \code{X}.
}
\value{
  A \code{function} with the same arguments as \code{f},
  which also belongs to the class \code{"funxy"}.
  This class has methods for
  \code{print}, \code{plot}, \code{contour} and \code{persp}.
}
\seealso{
  \code{\link{plot.funxy}},
  \code{\link{summary.funxy}}
}
\examples{
   f <- function(x,y) { x^2 + y^2 - 1} 
   g <- funxy(f, square(2))
   g
   ## evaluate function at any x, y coordinates
   g(0.2, 0.3)
   ## evaluate function at the points of a point pattern
   g(cells[1:4])
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{math}
