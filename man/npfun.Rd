\name{npfun}
\alias{npfun}
\title{
  Dummy Function Returns Number of Points
}
\description{
  Returns a summary function which is constant with value equal to
  the number of points in the point pattern.
}
\usage{
npfun(X, ..., r)
}
\arguments{
  \item{X}{
    Point pattern.
  }
  \item{\dots}{
    Ignored.
  }
  \item{r}{
    Vector of values of the distance argument \eqn{r}.
  }
}
\details{
  This function is normally not called by the user. Instead it is passed
  as an argument to the function \code{\link{psst}}.
}
\value{
  Object of class \code{"fv"} representing a constant function.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Ege Rubak and Jesper Moller.
}
\seealso{
  \code{\link{psst}}
}
\examples{
  data(cells)
  fit0 <- ppm(cells, ~1)
  v <- psst(fit0, npfun)
}
\keyword{spatial}
\keyword{nonparametric}

