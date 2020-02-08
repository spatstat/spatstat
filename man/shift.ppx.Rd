\name{shift.ppx}
\alias{shift.ppx}
\alias{shift.boxx}
\title{Apply Vector Translation To Box Or Point Pattern In Arbitrary Dimension}
\description{
  Applies a vector shift to a box or point pattern in arbitrary dimension
  (object of class \code{"boxx""} or \code{"ppx"}).
}
\usage{
 \method{shift}{boxx}(X, vec= 0, \dots)
 \method{shift}{ppx}(X, vec = 0, \dots, spatial = TRUE, temporal = TRUE, local = TRUE)
}
\arguments{
  \item{X}{Box or point pattern in arbitrary dimension (object of class \code{"boxx"} or \code{"ppx"}).}
  \item{vec}{Either a single numeric or a vector of the same length as the 
  dimension of the spatial and/or temporal and/or local domain.}
  \item{\dots}{Ignored}
  \item{spatial,temporal,local}{
    Logical to indicate whether or not to shift this type of coordinates for the
    \code{ppx} method.
  }
}
\value{
  For \code{shift.boxx}, another \code{"boxx"} object and for \code{shift.ppx}
  another \code{"ppx"} object. In both cases the new object represents the
  result of applying the vector shift.
}
\details{
  This is a method for the generic function \code{\link{shift}}.
}
\seealso{
  \code{\link{shift}},
  \code{\link{boxx}},
  \code{\link{ppx}}
}
\examples{
  vec <- c(2,3)
  dom <- boxx(c(0,1), c(0,1))
  X <- ppx(coords(cells), domain = dom)
  shift(dom, vec)
  Xs <- shift(X, vec)
  Xs
  head(coords(X), n = 3)
  head(coords(Xs), n = 3)
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{manip}
