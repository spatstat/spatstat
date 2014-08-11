\name{as.box3}
\Rdversion{1.1}
\alias{as.box3}
\title{
  Convert Data to Three-Dimensional Box
}
\description{
  Interprets data as the dimensions of a three-dimensional box.
}
\usage{
as.box3(...)
}
\arguments{
  \item{\dots}{
    Data that can be interpreted as giving the dimensions of a
    three-dimensional box. See Details.
  }
}
\details{
  This function converts data in various formats to
  an object of class \code{"box3"} representing a three-dimensional
  box (see \code{\link{box3}}). The arguments \code{\dots} may be
  \itemize{
    \item an object of class \code{"box3"}
    \item arguments acceptable to \code{box3}
    \item a numeric vector of length 6, interpreted as
    \code{c(xrange[1],xrange[2],yrange[1],yrange[2],zrange[1],zrange[2])}
    \item an object of class \code{"pp3"} representing a
    three-dimensional point pattern contained in a box.
  }
}
\value{
  Object of class \code{"box3"}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{box3}}, 
  \code{\link{pp3}}
}
\examples{
    X <- c(0,10,0,10,0,5)
    as.box3(X)
    X <- pp3(runif(42),runif(42),runif(42), box3(c(0,1)))
    as.box3(X)
}
\keyword{spatial}
\keyword{manip}
