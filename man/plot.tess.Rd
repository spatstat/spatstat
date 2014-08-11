\name{plot.tess}
\alias{plot.tess}
\title{Plot a tessellation}
\description{
  Plots a tessellation.
}
\usage{
  \method{plot}{tess}(x, ..., main, add=FALSE, col=NULL)
}
\arguments{
  \item{x}{Tessellation (object of class \code{"tess"}) to be plotted.}
  \item{\dots}{Arguments controlling the appearance of the plot.}
  \item{main}{Heading for the plot. A character string.}
  \item{add}{Logical. Determines whether the tessellation plot is added
    to the existing plot.
  }
  \item{col}{
    Colour of the tile boundaries. A character string.
    Ignored for pixel tessellations.
  }
    
}
\details{
  This is a method for the generic \code{\link{plot}} function
  for the class \code{"tess"} of tessellations (see \code{\link{tess}}).

  The arguments \code{\dots} control the appearance of the plot.
  They are passed to \code{\link{segments}},
  \code{\link{plot.owin}} or \code{\link{plot.im}}, depending on the
  type of tessellation. 
}
\value{
  None.
}
\seealso{
  \code{\link{tess}}
}
\examples{
  A <- tess(xgrid=0:4,ygrid=0:4)
  plot(A, col="blue", lwd=2, lty=2)
  B <- A[c(1, 2, 5, 7, 9)]
  plot(B, hatch=TRUE)
  v <- as.im(function(x,y){factor(round(5 * (x^2 + y^2)))}, W=owin())
  levels(v) <- letters[seq(length(levels(v)))]
  E <- tess(image=v)
  plot(E)
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{hplot}
