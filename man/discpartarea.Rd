\name{discpartarea}
\Rdversion{1.1}
\alias{discpartarea}
\title{
  Area of Part of Disc
}
\description{
  Compute area of intersection between a disc and a window  
}
\usage{
discpartarea(X, r, W=as.owin(X))
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"})
    specifying the centres of the discs.
    Alternatively, \code{X} may be in any format
    acceptable to \code{\link{as.ppp}}.
  }
  \item{r}{
    Matrix, vector or numeric value specifying the
    radii of the discs.
  }
  \item{W}{
    Window (object of class \code{"owin"}) with which the
    discs should be intersected.
  }
}
\details{
  This algorithm computes the exact area of the intersection between
  a window \code{W} and a disc (or each of several discs).
  The centres of the discs are specified by the point pattern
  \code{X}, and their radii are specified by \code{r}.

  If \code{r} is a single numeric value, then the algorithm computes the
  area of intersection between \code{W} and the disc of radius \code{r} centred
  at each point of \code{X}, and returns a one-column matrix
  containing one entry for each point of \code{X}.

  If \code{r} is a vector of length \code{m}, then the algorithm
  returns an \code{n * m} matrix in which the entry on row \code{i},
  column \code{j} is the area of the
  intersection between \code{W} and the disc centred at \code{X[i]}
  with radius \code{r[j]}.
  
  If \code{r} is a matrix, it should have one row for each point in
  \code{X}. The algorithm
  returns a matrix in which the entry on row \code{i},
  column \code{j} is the area of the
  intersection between \code{W} and the disc centred at \code{X[i]}
  with radius \code{r[i,j]}.

  Areas are computed by analytic geometry.
}
\value{
  Numeric matrix, with one row for each point of \code{X}.
}
\seealso{
  \code{\link{owin}},
  \code{\link{disc}}
}
\examples{
  data(letterR)
  X <- runifpoint(3, letterR)
  discpartarea(X, 0.2)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{math}
