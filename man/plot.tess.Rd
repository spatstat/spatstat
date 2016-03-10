\name{plot.tess}
\alias{plot.tess}
\title{Plot a tessellation}
\description{
  Plots a tessellation.
}
\usage{
  \method{plot}{tess}(x, ..., main, add=FALSE,
                      show.all=!add, col=NULL, do.plot=TRUE,
                      do.labels=FALSE,
                      labels=tilenames(x), labelargs=list())
}
\arguments{
  \item{x}{Tessellation (object of class \code{"tess"}) to be plotted.}
  \item{\dots}{Arguments controlling the appearance of the plot.}
  \item{main}{Heading for the plot. A character string.}
  \item{add}{Logical. Determines whether the tessellation plot is added
    to the existing plot.
  }
  \item{show.all}{
    Logical value indicating whether to plot everything
    including the main title and the observation window of \code{x}.
  }
  \item{col}{
    Colour of the tile boundaries. A character string.
    Ignored for pixel tessellations.
  }
  \item{do.plot}{
    Logical value indicating whether to actually perform the plot.
  }
  \item{do.labels}{
    Logical value indicating whether to show a text label for each tile
    of the tessellation.
  }
  \item{labels}{Character vector of labels for the tiles.}
  \item{labelargs}{
    List of arguments passed to
    \code{\link{text.default}} to control display of the text labels.
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
  (Invisible) window of class \code{"owin"} specifying a bounding box
  for the plot (including a colour ribbon if plotted).
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

\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{hplot}
