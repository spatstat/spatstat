\name{convexhull}
\alias{convexhull}
\title{Convex Hull}
\description{
  Computes the convex hull of a spatial object.
}
\usage{
convexhull(x)
}
\arguments{
  \item{x}{
    a window (object of class \code{"owin"}),
    a point pattern (object of class \code{"ppp"}),
    a line segment pattern (object of class \code{"psp"}),
    or an object that can be converted to a window
    by \code{\link{as.owin}}.
  }
}
\value{
  A window (an object of class \code{"owin"}).
}
\details{
  This function computes the convex hull of the spatial object \code{x}.
}
\seealso{
  \code{\link{owin}},
  \code{\link{convexhull.xy}},
  \code{\link{is.convex}}
}
\examples{
   data(demopat)
   W <- Window(demopat)
   plot(convexhull(W), col="lightblue", border=NA)
   plot(W, add=TRUE, lwd=2)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{utilities}




