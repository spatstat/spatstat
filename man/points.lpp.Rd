\name{points.lpp}
\alias{points.lpp}
\title{
  Draw Points on Existing Plot
}
\description{
  For a point pattern on a linear network, this function draws the
  coordinates of the points only, on the existing plot display.
}
\usage{
\method{points}{lpp}(x, \dots)
}
\arguments{
  \item{x}{
    A point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{\dots}{
    Additional arguments passed to \code{\link[graphics]{points.default}}.
  }
}
\details{
  This is a method for the generic function
  \code{\link[graphics]{points}}
  for the class \code{"lpp"} of point patterns on a linear network.

  If \code{x} is a point pattern on a linear network, then
  \code{points(x)} plots the spatial coordinates of the points only,
  on the existing plot display, 
  without plotting the underlying network.
  It is an error to call this function if a plot has not yet been
  initialised.
  
  The spatial coordinates are extracted and passed to
  \code{\link[graphics]{points.default}} along with any extra arguments.
  Arguments controlling the colours and the plot symbols are interpreted
  by \code{\link[graphics]{points.default}}. For example, if the
  argument \code{col} is a vector, then the \code{i}th point is drawn
  in the colour \code{col[i]}.
}
\section{Difference from plot method}{
  The more usual way to plot the points is using \code{\link{plot.lpp}}.
  For example \code{plot(x)} would plot both the points and the
  underlying network, while \code{plot(x, add=TRUE)} would plot only the
  points. The interpretation of arguments controlling the colours and
  plot symbols is different here: they determine a symbol map, as explained
  in the help for \code{\link{plot.ppp}}.
}
\value{
  Null.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{plot.lpp}}, \code{\link[graphics]{points.default}}
}
\examples{
 plot(Frame(spiders), main="Spiders on a Brick Wall")
 points(spiders)
}
\keyword{spatial}
\keyword{hplot}

