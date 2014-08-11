\name{areaGain}
\alias{areaGain}
\title{Difference of Disc Areas}
\description{
  Computes the area of that part of a disc
  that is not covered by other discs.
}
\usage{
   areaGain(u, X, r, ..., W=as.owin(X), exact=FALSE,
                     ngrid=spatstat.options("ngrid.disc"))
}
\arguments{
  \item{u}{
    Coordinates of the centre of the disc of interest.
    A vector of length 2.
    Alternatively, a point pattern (object of class \code{"ppp"}).
  }
  \item{X}{
    Locations of the centres of other discs.
    A point pattern (object of class \code{"ppp"}).
  }
  \item{r}{
    Disc radius, or vector of disc radii.
  }
  \item{\dots}{Ignored.}
  \item{W}{
    Window (object of class \code{"owin"})
    in which the area should be computed.
  }
  \item{exact}{
    Choice of algorithm.
    If \code{exact=TRUE}, areas are computed exactly using
    analytic geometry. If \code{exact=FALSE} then a faster algorithm
    is used to compute a discrete approximation to the areas.
  }
  \item{ngrid}{
    Integer. Number of points in the square grid used to compute
    the discrete approximation, when \code{exact=FALSE}.
}
}
\value{
  A matrix with one row for each point in \code{u}
  and one column for each value in \code{r}. 
}
\details{
  This function computes the area of that part of
  the disc of radius \code{r} centred at the location \code{u}
  that is \emph{not} covered by any of the discs of radius \code{r}
  centred at the points of the pattern \code{X}.
  This area is important in some calculations related to
  the area-interaction model \code{\link{AreaInter}}.

  If \code{u} is a point pattern and \code{r} is a vector,
  the result is a matrix, with one row for each point in \code{u}
  and one column for each entry of \code{r}. The \code{[i,j]} entry
  in the matrix is the area of that part of the disc of radius
  \code{r[j]} centred at the location \code{u[i]} that is
  \emph{not} covered by any of the discs of radius \code{r[j]}
  centred at the points of the pattern \code{X}.

  If \code{W} is not \code{NULL}, then the areas are computed only
  inside the window \code{W}. 
}
\seealso{
  \code{\link{AreaInter}},
  \code{\link{areaLoss}}
}
\examples{
   data(cells)
   u <- c(0.5,0.5)
   areaGain(u, cells, 0.1)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
