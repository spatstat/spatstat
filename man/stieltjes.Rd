\name{stieltjes}
\alias{stieltjes}
\title{Compute Integral of Function Against Cumulative Distribution}
\description{
  Computes the Stieltjes integral 
  of a function \eqn{f} with respect to a function \eqn{M}.
}
\usage{
stieltjes(f, M, ...)
}
\arguments{
  \item{f}{
    The integrand. A function in the \R language.
  }
  \item{M}{
    The cumulative function against which \code{f} will be
    integrated. An object of class \code{"fv"}.
  }
  \item{\dots}{
    Additional arguments passed to \code{f}.
  }
}
\details{
  This command computes the Stieltjes integral
  \deqn{I = \int f(x) dM(x)}{I = integral f(x) dM(x)}
  of a real-valued function \eqn{f(x)}
  with respect to a nondecreasing function \eqn{M(x)}.

  One common use of the Stieltjes integral is
  to find the mean value of a random variable from its
  cumulative distribution function \eqn{F(x)}. The mean value is
  the Stieltjes integral of \eqn{f(x)=x} with respect to \eqn{F(x)}.

  The argument \code{f} should be a \code{function} in the \R language.
  It should accept a numeric vector argument \code{x} and should return
  a numeric vector of the same length.

  The argument \code{M} should be a function value table
  (object of class \code{"fv"}, see \code{\link{fv.object}}).
  Such objects are returned
  by the commands \code{\link{Kest}}, \code{\link{Gest}}, etc.
}
\value{
  A list containing the value of the Stieltjes integral
  computed using each of the versions of the function \code{M}.
}
\seealso{
  \code{\link{fv.object}},
  \code{\link{Gest}}
}
\examples{
  data(redwood)
  # estimate cdf of nearest neighbour distance
  G <- Gest(redwood)
  # compute estimate of mean nearest neighbour distance
  stieltjes(function(x){x}, G)
  # estimated probability of a distance in the interval [0.1,0.2]
  stieltjes(function(x,a,b){ (x >= a) & (x <= b)}, G, a=0.1, b=0.2)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
