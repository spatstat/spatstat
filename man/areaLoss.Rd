\name{areaLoss}
\alias{areaLoss}
\title{Difference of Disc Areas}
\description{
  Computes the area of that part of a disc
  that is not covered by other discs.
}
\usage{
   areaLoss(X, r, ..., W=as.owin(X), subset=NULL,
                 exact=FALSE,
                 ngrid=spatstat.options("ngrid.disc"))
}
\arguments{
  \item{X}{
    Locations of the centres of discs.
    A point pattern (object of class \code{"ppp"}).
  }
  \item{r}{
    Disc radius, or vector of disc radii.
  }
  \item{\dots}{Ignored.}
  \item{W}{
    Optional. Window (object of class \code{"owin"})
    inside which the area should be calculated.
  }
  \item{subset}{
    Optional. Index identifying a subset of the points of \code{X}
    for which the area difference should be computed.
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
  A matrix with one row for each point in \code{X} (or \code{X[subset]})
  and one column for each value in \code{r}. 
}
\details{
  This function computes, for each point \code{X[i]} in \code{X}
  and for each radius \code{r},
  the area of that part of the disc of radius \code{r} centred at the
  location \code{X[i]} that is \emph{not} covered by any of the
  other discs of radius \code{r} centred at the points \code{X[j]}
  for \code{j} not equal to \code{i}.
  This area is important in some calculations related to
  the area-interaction model \code{\link{AreaInter}}.

  The result is a matrix, with one row for each point in \code{X}
  and one column for each entry of \code{r}. 
}
\seealso{
  \code{\link{AreaInter}},
  \code{\link{areaGain}},
  \code{\link{dilated.areas}}
}
\examples{
   data(cells)
   areaLoss(cells, 0.1)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{math}
