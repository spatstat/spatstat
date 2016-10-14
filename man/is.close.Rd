\name{is.close}
\alias{is.close}
\alias{is.close.default}
\alias{is.close.ppp}
\alias{is.close.pp3}
\title{
  Check Whether Points Have Close Neighbours
}
\description{
  For each point in a point pattern, determine whether the
  point has a close neighbour in the same pattern.
}
\usage{
  is.close(X, r, Y=NULL)

  \method{is.close}{default}(X,r, Y=NULL)

  \method{is.close}{ppp}(X,r, Y=NULL)

  \method{is.close}{pp3}(X,r, Y=NULL)
}
\arguments{
  \item{X,Y}{
    Point patterns of class \code{"ppp"} or \code{"pp3"} or \code{"lpp"}.
  }
  \item{r}{
    Threshold distance: a number greater than zero.
  }
}
\details{
  This is simply a faster version of \code{(nndist(X) <= r)}
  or \code{(nncross(X,Y,what="dist") <= r)}.

  \code{is.close(X,r)} determines, for each point in the pattern \code{X},
  whether or not this point has a neighbour in the same pattern \code{X}
  which lies at a distance less than or equal to \code{r}.

  \code{is.close(X,r,Y)} determines, for each point in the pattern \code{X},
  whether or not this point has a neighbour in the \emph{other} pattern
  \code{Y} which lies at a distance less than or equal to \code{r}.

  The function \code{is.close} is generic, with methods for
  \code{"ppp"} and \code{"pp3"} and a default method.
}
\value{
  A logical vector, with one entry for each point of \code{X}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{nndist}}
}
\examples{
  is.close(redwood, 0.05)
  with(split(amacrine), is.close(on, 0.05, off))
}
\keyword{spatial}
\keyword{math}

