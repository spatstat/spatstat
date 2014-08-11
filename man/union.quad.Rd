\name{union.quad}
\alias{union.quad}
\title{Union of Data and Dummy Points}
\description{
  Combines the data and dummy points of a quadrature scheme
  into a single point pattern.
}
\usage{
 union.quad(Q)
}
\arguments{
  \item{Q}{A quadrature scheme (an object of class \code{"quad"}).}
}
\value{
  A point pattern (of class \code{"ppp"}).
}
\details{
  The argument \code{Q} should be a quadrature scheme (an object of class
  \code{"quad"}, see \code{\link{quad.object}} for details).
  
  This function combines the data and dummy points of \code{Q}
  into a single point pattern. If either the data or the dummy points
  are marked, the result is a marked point pattern.

  The function \code{\link{as.ppp}} will perform the same task.
}
\seealso{
  \code{\link{quad.object}},
  \code{\link{as.ppp}}
}
\examples{
  data(simdat)
  Q <- quadscheme(simdat, default.dummy(simdat))
  U <- union.quad(Q)
  \dontrun{plot(U)}
  # equivalent:
  U <- as.ppp(Q)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
