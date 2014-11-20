\name{yardstick}
\alias{yardstick}
\title{
  Yardstick or Scale Bar 
}
\description{
  Create a spatial object that represents a scale bar.
}
\usage{
yardstick(x0, y0, x1, y1, txt)
}
\arguments{
  \item{x0,y0,x1,y1}{
    Spatial coordinates of both ends of the line segment.
    Alternatively \code{x0} can be a point pattern (class \code{"ppp"})
    containing exactly two points, or a line segment pattern
    (class \code{"psp"}) consisting of exactly one line segment.
  }
  \item{txt}{
    The text to be displayed beside the line segment.
    Either a character string or an expression.
  }
}
\details{
  A yardstick or scale bar is a line segment, drawn on any spatial
  graphics display, indicating the scale of the plot. 

  This command creates an object of class \code{"yardstick"}
  which is essentially a line segment joining the
  endpoints \code{(x0,y0)} and \code{(x1,y1)} and annotated by
  the text \code{txt}. Plotting this object will display the
  yardstick on the current plot, with the annotation next to the line.

  A yardstick can be one of the components in a
  list of spatial objects (see \code{\link{solist}})
  or in a layered object (see \code{\link{layered}}).
}
\value{
  Object of class \code{"yardstick"} which can be plotted.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\seealso{
  \code{\link{plot.yardstick}}
}
\examples{
  X <- rescale(swedishpines)
  plot(X, pch=16, main="")
  ys <- yardstick(as.psp(list(xmid=4, ymid=0.5, length=1, angle=0),
                         window=Window(X)),
                  txt="1 m")
  plot(ys)
}
\keyword{spatial}
\keyword{hplot}

