\name{pppdist}
\alias{pppdist}
\title{Distance Between Two Point Patterns}
\description{
  Given two point patterns, find the distance between them based on
  optimal point matching.
}
\usage{
  pppdist(X, Y, type = "spa", cutoff = 1, q = 1, matching = TRUE,
    ccode = TRUE, precision = NULL, approximation = 10,
    show.rprimal = FALSE, timelag = 0)
}
\arguments{
  \item{X,Y}{Two point patterns (objects of class \code{"ppp"}).}
  \item{type}{
    A character string giving the type of distance to be computed.
    One of \code{"spa"} (default), \code{"ace"} or \code{"mat"}, indicating
    whether the algorithm should find the optimal matching based on "subpattern
    assignment", "assignment only if cardinalities are equal" or "mass transfer".
    See details below. 
  }
  \item{cutoff}{
    The value \eqn{> 0} at which interpoint distances are cut off.
  }
  \item{q}{
    The order of the average that is applied to the interpoint distances.
    May be \code{Inf}, in which case the maximum of the interpoint distances is taken.
  }
  \item{matching}{
    Logical. Whether to return the optimal matching or only the associated distance.
  }
  \item{ccode}{
    Logical. If \code{FALSE}, \R code is used which allows for higher precision, but is
    much slower.
  }
  \item{precision}{
    Index controlling accuracy of algorithm. The \code{q}-th powers of interpoint distances
    will be rounded to the nearest multiple of \code{10^(-precision)}. There is a sensible default which depends on \code{ccode}.
  }
  \item{approximation}{
    If \code{q = Inf}, compute distance based on the optimal matching for the
    corresponding distance of order \code{approximation}. Can be \code{Inf}, but
    this makes computations extremely slow.
  }
  \item{show.rprimal}{
    Logical. Whether to display a plot showing the
    iterative solution of the restricted primal problem.
  }
  \item{timelag}{
    Time lag, in seconds, between successive displays of the
    iterative solution of the restricted primal problem.
  }
}
\details{
  Computes the distance between point patterns \code{X} and \code{Y} based
  on finding the matching between them
  which minimizes the average of the distances between matched points
  (if \code{q=1}), the maximum distance between matched points
  (if \code{q=Inf}), and in general the \code{q}-th order average
  (i.e. the \code{1/q}th power of the sum of
  the \code{q}th powers) of the distances between matched points.
  Distances between matched points are Euclidean distances cut off at
  the value of \code{cutoff}.

  The parameter \code{type} controls the behaviour of the algorithm if
  the cardinalities of the point patterns are different. For the type \code{"spa"}
  (subpattern assignment) the subpattern of the point pattern
  with the larger cardinality \eqn{n} that is closest to the point pattern
  with the smaller cardinality \eqn{m} is determined; then the \code{q}-th order
  average is taken over \eqn{n} values: the \eqn{m} distances of matched points
  and \eqn{n-m} "penalty distances" of value \code{cutoff} for
  the unmatched points. For the type \code{"ace"} (assignment only if 
  cardinalities equal) the matching is empty and the distance returned is equal
  to \code{cutoff} if the cardinalities differ. For the
  type \code{"mat"} (mass transfer) each point pattern is assumed
  to have total mass \eqn{m} (= the smaller cardinality) distributed evenly
  among its points; the algorithm finds then the "mass transfer plan" that
  minimizes the \code{q}-th order weighted average of the distances, where 
  the weights are given by the transferred mass divided by \eqn{m}. The result is a fractional matching (each match of two points has a weight in \eqn{(0,1]}) 
  with the minimized quantity as the associated distance.
  
  The computations for all three types rely heavily on a specialized
  primal-dual algorithm (described in Luenberger (2003), Section 5.9)
  for Hitchcock's problem of optimal transport of a product from a number
  of suppliers to a number of (e.g. vending) locations. The C implementation
  used by default can handle patterns with a few hundreds of points, but
  should not be used with thousands of points. By setting \code{show.rprimal = TRUE},
  some insight in the working of the algorithm can be gained.
  
  For moderate and large values of \code{q} there
  can be numerical issues based on the fact that the \code{q}-th powers of
  distances are taken and some positive values enter the optimization algorithm
  as zeroes because they are too small in comparison with the larger values.
  In this case the number of zeroes introduced is given in a warning message,
  and it is possible then that the matching obtained is not optimal and the associated
  distance is only a strict upper bound of the true distance.
  As a general guideline (which can be very wrong in special situations) a small
  number of zeroes (up to about 50\% of the smaller point pattern
  cardinality \eqn{m})
  usually still results in the right matching, and the number can even be quite a bit
  higher and usually still provides a highly accurate upper bound for the distance.
  These numerical
  problems can be reduced by enforcing (much slower) \R code via the 
  argument \code{ccode = FALSE}. 

  For \code{q = Inf} there is no fast algorithm available, which is why approximation is
  normally used: for finding the optimal matching, \code{q} is
  set to the value of \code{approximation}. The
  resulting distance is still given as the maximum rather than the
  \code{q}-th order average in the corresponding distance computation.
  If \code{approximation = Inf}, approximation is suppressed and a very inefficient
  exhaustive search for the best matching is performed.

  The value of \code{precision} should normally not be supplied by the user. If
  \code{ccode = TRUE}, this value is preset to the highest exponent of 10 that
  the C code still can handle (usually \eqn{9}). If \code{ccode = FALSE}, the value is
  preset according to \code{q} (usually \eqn{15} if \code{q} is small),
  which can sometimes be changed to obtain less severe warning messages.
}
\value{
  Normally an object of class \code{pppmatching} that contains detailed
  information about the parameters used and the resulting distance.
  See \code{\link{pppmatching.object}} for details.
  If \code{matching = FALSE}, only the numerical value of the distance
  is returned.
}
\references{
  Hitchcock F.L. (1941)
  The distribution of a product from several sources to numerous
  localities. \emph{J. Math. Physics} \bold{20}, 224--230.

  Luenberger D.G. (2003). \emph{Linear and nonlinear programming.}
  Second edition. Kluwer.

  Schuhmacher, D. and Xia, A. (2008)
  A new metric between distributions of point processes.
  \emph{Advances in Applied Probability} \bold{40}, 651--672

  Schuhmacher, D., Vo, B.-T. and Vo, B.-N. (2008)
  A consistent metric for performance evaluation of multi-object
  filters.
  \emph{IEEE Transactions on Signal Processing} \bold{56}, 3447--3457.
}
\author{
  Dominic Schuhmacher
  \email{dominic.schuhmacher@stat.unibe.ch}
  \url{http://www.dominic.schuhmacher.name}
}
\seealso{
  \code{\link{pppmatching.object}}
  \code{\link{matchingdist}}
}
\examples{
  # equal cardinalities
  X <- runifpoint(100)
  Y <- runifpoint(100)
  m <- pppdist(X, Y)
  m
  \dontrun{
  plot(m)
  }
  
  # differing cardinalities
  X <- runifpoint(14)
  Y <- runifpoint(10)
  m1 <- pppdist(X, Y, type="spa")
  m2 <- pppdist(X, Y, type="ace")
  m3 <- pppdist(X, Y, type="mat")
  summary(m1)
  summary(m2)
  summary(m3)
  \dontrun{
  m1$matrix
  m2$matrix
  m3$matrix
  }

  # q = Inf
  X <- runifpoint(10)
  Y <- runifpoint(10)
  mx1 <- pppdist(X, Y, q=Inf)$matrix
  mx2 <- pppdist(X, Y, q=Inf, ccode=FALSE, approximation=50)$matrix
  mx3 <- pppdist(X, Y, q=Inf, approximation=Inf)$matrix
  ((mx1 == mx2) && (mx2 == mx3))
       # TRUE if approximations are good
}
\keyword{spatial}
\keyword{math}
