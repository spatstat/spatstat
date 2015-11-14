\name{overlap.owin}
\alias{overlap.owin}
\title{
  Compute Area of Overlap
}
\description{
  Computes the area of the overlap (intersection) of
  two windows.
}
\usage{
overlap.owin(A, B)
}
\arguments{
  \item{A,B}{
    Windows (objects of class \code{"owin"}).
  }
}
\details{
  This function computes the area of the overlap between the two
  windows \code{A} and \code{B}.

  If one of the windows is a binary mask, then both windows
  are converted to masks on the same grid, and the area is
  computed by counting pixels. Otherwise, the area is computed
  analytically (using the discrete Stokes theorem).
}
\value{
  A single numeric value.
}
\seealso{
  \code{\link{intersect.owin}}, 
  \code{\link{area.owin}},
  \code{\link{setcov}}.
}
\examples{
  A <- square(1)
  B <- shift(A, c(0.3, 0.2))
  overlap.owin(A, B)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{math}