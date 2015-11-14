\name{yardstick}
\alias{textstring}
\alias{onearrow}
\alias{yardstick}
\title{
  Text, Arrow or Scale Bar in a Diagram
}
\description{
  Create spatial objects that represent a text string,
  an arrow, or a yardstick (scale bar).
}
\usage{
textstring(x, y, txt = NULL, \dots)

onearrow(x0, y0, x1, y1, txt = NULL, \dots)

yardstick(x0, y0, x1, y1, txt = NULL, \dots)
}
\arguments{
  \item{x,y}{
    Coordinates where the text should be placed.
  }
  \item{x0,y0,x1,y1}{
    Spatial coordinates of both ends of the arrow or yardstick.
    Alternatively \code{x0} can be a point pattern (class \code{"ppp"})
    containing exactly two points, or a line segment pattern
    (class \code{"psp"}) consisting of exactly one line segment.
  }
  \item{txt}{
    The text to be displayed beside the line segment.
    Either a character string or an expression.
  }
  \item{\dots}{
    Additional named arguments for plotting the object.
  }
}
\details{
  These commands create objects that represent components of a
  diagram:
  \itemize{
    \item \code{textstring} creates an object that represents a string
    of text at a particular spatial location.
    \item \code{onearrow} creates an object that represents an arrow
    between two locations.
    \item \code{yardstick} creates an object that represents
    a scale bar: a line segment indicating the scale of the plot.
  }
  To display the relevant object, it should be plotted, using
  \code{plot}. See the help files for the plot methods
  \code{\link{plot.textstring}}, \code{\link{plot.onearrow}}
  and \code{\link{plot.yardstick}}.

  These objects are designed to be included as components in a
  \code{\link{layered}} object or a \code{\link{solist}}. This makes it
  possible to build up a diagram consisting of many spatial objects,
  and to annotate the diagram with arrows, text and so on, so that
  ultimately the entire diagram is plotted using \code{plot}.
}
\value{
  An object of class \code{"diagramobj"} which also
  belongs to one of the special classes \code{"textstring"}, \code{"onearrow"}
  or \code{"yardstick"}. There are methods for \code{plot},
  \code{print}, \code{"["} and \code{\link{shift}}.
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\seealso{
  \code{\link{plot.textstring}},
  \code{\link{plot.onearrow}},
  \code{\link{plot.yardstick}}.
}
\examples{
  X <- rescale(swedishpines)
  plot(X, pch=16, main="")
  ys <- yardstick(as.psp(list(xmid=4, ymid=0.5, length=1, angle=0),
                         window=Window(X)),
                  txt="1 m")
  plot(ys, angle=90)
}
\keyword{spatial}
\keyword{hplot}

