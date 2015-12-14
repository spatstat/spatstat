\name{distfun.lpp}
\Rdversion{1.1}
\alias{distfun.lpp}
\title{
  Distance Map on Linear Network
}
\description{
  Compute the distance function of a point pattern on a linear network.
}
\usage{
  \method{distfun}{lpp}(X, ..., k=1)
}
\arguments{
  \item{X}{
    A point pattern on a linear network
    (object of class \code{"lpp"}).
  }
  \item{k}{
    An integer. The distance to the \code{k}th nearest point
    will be computed.
  }
  \item{\dots}{
    Extra arguments are ignored. 
  }
}
\details{
  On a linear network \eqn{L}, the \dQuote{geodesic distance function}
  of a set of points \eqn{A} in \eqn{L} is the
  mathematical function \eqn{f} such that, for any 
  location \eqn{s} on \eqn{L},
  the function value \code{f(s)}
  is the shortest-path distance from \eqn{s} to \eqn{A}.

  The command \code{distfun.lpp} is a method for the generic command
  \code{\link{distfun}}
  for the class \code{"lpp"} of point patterns on a linear network.

  If \code{X} is a point pattern on a linear network,
  \code{f <- distfun(X)} returns a \emph{function}
  in the \R language that represents the
  distance function of \code{X}. Evaluating the function \code{f}
  in the form \code{v <- f(x,y)}, where \code{x} and \code{y}
  are any numeric vectors of equal length containing coordinates of
  spatial locations, yields the values of the distance function at these
  locations. More efficiently \code{f} can take the arguments
  \code{x, y, seg, tp} where \code{seg} and \code{tp} are the local
  coordinates on the network.

  The function \code{f} obtained from \code{f <- distfun(X)}
  also belongs to the class \code{"linfun"}. 
  It can be printed and plotted immediately as shown in the Examples.
  It can be 
  converted to a pixel image using \code{\link{as.linim}}. 
}
\value{
  A \code{function} with arguments \code{x,y} and optional
  arguments \code{seg,tp}.
  It also belongs to the class \code{"linfun"} which has methods
  for \code{plot}, \code{print} etc.
}
\seealso{
  \code{\link{linfun}},
  \code{\link{methods.linfun}}.

  To identify \emph{which} point is the nearest neighbour, see
  \code{\link{nnfun.lpp}}.
}
\examples{
   data(letterR)
   X <- runiflpp(3, simplenet)
   f <- distfun(X)
   f
   plot(f)

   # using a distfun as a covariate in a point process model:
   Y <- runiflpp(4, simplenet)
   fit <- lppm(Y, ~D, covariates=list(D=f))
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{math}
