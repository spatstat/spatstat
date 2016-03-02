\name{pppmatching}
\alias{pppmatching}
\title{Create a Point Matching}
\description{
  Creates an object of class \code{"pppmatching"} representing
  a matching of two planar point patterns (objects of class \code{"ppp"}).
}
\usage{
  pppmatching(X, Y, am, type = NULL, cutoff = NULL, q = NULL,
    mdist = NULL)
}
\arguments{
  \item{X,Y}{Two point patterns (objects of class \code{"ppp"}).}
  \item{am}{
    An \code{npoints(X)} by \code{npoints(Y)} matrix with entries \eqn{\geq 0}{>= 0} 
    that specifies which points are matched and with what weight;
    alternatively, an object that can be coerced to this form
    by \code{as.matrix}.
  }
  \item{type}{
    A character string giving the type of the matching.
    One of \code{"spa"}, \code{"ace"} or \code{"mat"}, or \code{NULL}
    for a generic or unknown matching.
  }
  \item{cutoff, q}{
    Numerical values specifying the cutoff value \eqn{> 0} for interpoint distances and
    the order \eqn{q \in [1,\infty]}{q in [0,Inf]} of the average that is applied to them.
    \code{NULL} if not applicable or unknown.
  }
  \item{mdist}{
    Numerical value for the distance to be associated with the matching.
  }
}
\details{
  The argument \code{am} is interpreted as a "generalized adjacency matrix":
  if the \code{[i,j]}-th entry is positive, then the \code{i}-th point
  of \code{X} and the \code{j}-th point of \code{Y} are matched and the
  value of the entry gives the corresponding weight of the match. For
  an unweighted matching all the weights should be set to \eqn{1}.

  The remaining arguments are optional and allow to save
  additional information about the matching. See the help files for
  \code{\link{pppdist}} and \code{\link{matchingdist}} for details on
  the meaning of these parameters.
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
  # a random unweighted complete matching
  X <- runifpoint(10)
  Y <- runifpoint(10)
  am <- r2dtable(1, rep(1,10), rep(1,10))[[1]]
        # generates a random permutation matrix
  m <- pppmatching(X, Y, am)
  summary(m)
  m$matrix
  \dontrun{
    plot(m)
  }

  # a random weighted complete matching
  X <- runifpoint(7)
  Y <- runifpoint(7)
  am <- r2dtable(1, rep(10,7), rep(10,7))[[1]]/10
        # generates a random doubly stochastic matrix
  m2 <- pppmatching(X, Y, am)
  summary(m2)
  m2$matrix
  \dontrun{
    # Note: plotting does currently not distinguish
    # between different weights
    plot(m2)
  }
}
\keyword{spatial}
\keyword{datagen}
