\name{matchingdist}
\alias{matchingdist}
\title{Distance for a Point Pattern Matching}
\description{
  Computes the distance associated with a matching between two point patterns.
}
\usage{
  matchingdist(matching, type = NULL, cutoff = NULL, q = NULL)
}
\arguments{
  \item{matching}{A point pattern matching (an object of class \code{"pppmatching"}).}
  \item{type}{
    A character string giving the type of distance to be computed.
    One of \code{"spa"}, \code{"ace"} or \code{"mat"}. See details below. 
  }
  \item{cutoff}{
    The value \eqn{> 0} at which interpoint distances are cut off.
  }
  \item{q}{
    The order of the average that is applied to the interpoint distances.
    May be \code{Inf}, in which case the maximum of the interpoint distances is taken.
  }
}
\details{
  Computes the distance specified by \code{type}, \code{cutoff}, and \code{order}
  for a point matching. If any of these arguments are not provided, the function
  uses the corresponding elements of \code{matching} (if available).

  For the type \code{"spa"} (subpattern assignment) it is assumed that the points
  of the point pattern with the smaller cardinality \eqn{m} are matched to a
  \eqn{m}-point subpattern of the point pattern with the larger
  cardinality \eqn{n} in a 1-1 way. The distance
  is then given as the \code{q}-th order average of the \eqn{m} distances between
  matched points (minimum of Euclidean distance and \code{cutoff})
  and \eqn{n-m} "penalty distances" of value \code{cutoff}.

  For the type \code{"ace"} (assignment only if cardinalities equal) the matching
  is assumed to be 1-1 if the cardinalities of the point patterns are
  the same, in which case the \code{q}-th order average of the matching distances
  (minimum of Euclidean distance and \code{cutoff}) is taken. If the cardinalities
  are different, the matching may be arbitrary and the distance returned is always
  equal to \code{cutoff}.

  For the type \code{mat} (mass transfer) it is assumed that each point of
  the point pattern with the smaller cardinality \eqn{m} has mass \eqn{1},
  each point of the point pattern with the larger cardinality \eqn{n}
  has mass \eqn{m/n},
  and fractions of these masses are matched in such a way that each point
  contributes exactly its mass. The distance is then given as the \code{q}-th
  order weighted average of all distances (minimum of Euclidean distance
  and \code{cutoff}) of (partially) matched points with weights equal to the
  fractional masses divided by \eqn{m}.

  If the cardinalities of the two point patterns are equal,
  \code{matchingdist(m, type, cutoff, q)} yields the same result
  no matter if \code{type} is \code{"spa"}, \code{"ace"} or
  \code{"mat"}.
}
\value{
  Numeric value of the distance associated with the matching.
}
\author{
  Dominic Schuhmacher
  \email{dominic.schuhmacher@stat.unibe.ch}
  \url{http://www.dominic.schuhmacher.name}
}
\seealso{
  \code{\link{pppdist}}
  \code{\link{pppmatching.object}}
}
\examples{
  # an optimal matching
  X <- runifpoint(20)
  Y <- runifpoint(20)
  m.opt <- pppdist(X, Y)
  summary(m.opt)
  matchingdist(m.opt)
       # is the same as the distance given by summary(m.opt)
  
  # sequential nearest neighbour matching
  # (go through all points of point pattern X in sequence
  # and match each point with the closest point of Y that is
  # still unmatched)
  am <- matrix(0, 20, 20)
  h <- matrix(c(1:20, rep(0,20)), 20, 2)
  h[1,2] = nncross(X[1],Y)[1,2]
  for (i in 2:20) {
    nn <- nncross(X[i],Y[-h[1:(i-1),2]])[1,2]
    h[i,2] <- ((1:20)[-h[1:(i-1),2]])[nn]
  }
  am[h] <- 1
  m.nn <- pppmatching(X, Y, am)
  matchingdist(m.nn, type="spa", cutoff=1, q=1)
       # is >= the distance obtained for m.opt
       # in most cases strictly >

  opa <- par(mfrow=c(1,2))
  plot(m.opt, main="optimal")
  plot(m.nn,  main="nearest neighbour")
  text(X, 1:20, pos=1, offset=0.3, cex=0.8)
  par(opa)
}
\keyword{spatial}
\keyword{math}
