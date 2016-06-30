\name{as.solist}
\alias{as.solist}
\title{
  Convert List of Two-Dimensional Spatial Objects
}
\description{
  Given a list of two-dimensional spatial objects,
  convert it to the class \code{"solist"}.
}
\usage{
as.solist(x, \dots)
}
\arguments{
  \item{x}{
    A list of objects, each representing a two-dimensional
    spatial dataset.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{solist}}.
  }
}
\details{
  This command makes the list \code{x} into an object of class \code{"solist"}
  (spatial object list).
  See \code{\link{solist}} for details.

  The entries in the list \code{x} should be
  two-dimensional spatial datasets
  (not necessarily of the same class).
}
\value{
  A list, usually of class \code{"solist"}.
}
\seealso{
  \code{\link{solist}},
  \code{\link{as.anylist}},
  \code{\link{solapply}}.
}
\examples{
  x <- list(cells, density(cells))
  y <- as.solist(x)
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\keyword{spatial}
\keyword{list}
\keyword{manip}
