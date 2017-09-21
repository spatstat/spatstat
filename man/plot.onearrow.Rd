\name{plot.onearrow}
\alias{plot.onearrow}
\title{Plot an Arrow}
\description{Plots an object of class \code{"onearrow"}.}
\usage{
 \method{plot}{onearrow}(x, \dots,
  add = FALSE, main = "",
  retract = 0.05, headfraction = 0.25, headangle = 12, headnick = 0.1,
  col.head = NA, lwd.head = lwd, lwd = 1, col = 1,
  zap = FALSE, zapfraction = 0.07,
  pch = 1, cex = 1, do.plot = TRUE, do.points = FALSE, show.all = !add)
}
\arguments{
  \item{x}{
    Object of class \code{"onearrow"} to be plotted.
    This object is created by the command \code{\link{onearrow}}.
  }
  \item{\dots}{
    Additional graphics arguments passed to
    \code{\link[graphics]{segments}} to control the appearance of the line.
  }
  \item{add}{Logical value indicating whether to add graphics to the
    existing plot (\code{add=TRUE}) or to start a new plot
    (\code{add=FALSE}).
  }
  \item{main}{Main title for the plot.}
  \item{retract}{
    Fraction of length of arrow to remove at each end.
  }
  \item{headfraction}{
    Length of arrow head as a fraction of overall length of arrow.
  }
  \item{headangle}{
    Angle (in degrees) between the outer edge of the arrow head
    and the shaft of the arrow.
  }
  \item{headnick}{
    Size of the nick in the trailing edge of the arrow head
    as a fraction of length of arrow head. 
  }
  \item{col.head,lwd.head}{
    Colour and line style of the filled arrow head.
  }
  \item{col,lwd}{
    Colour and line style of the arrow shaft.
  }
  \item{zap}{
    Logical value indicating whether the arrow should
    include a Z-shaped (lightning-bolt) feature in the
    middle of the shaft.
  }
  \item{zapfraction}{
   Size of Z-shaped deviation as a fraction of total arrow length.
  }
  \item{pch,cex}{
    Plot character and character size for the two end points of the arrow,
    if \code{do.points=TRUE}.
  }
  \item{do.plot}{
    Logical. Whether to actually perform the plot.
  }
  \item{do.points}{
    Logical. Whether to display the two end points of the arrow as well.
  }
  \item{show.all}{
    Internal use only.
  }
}
\details{
  The argument \code{x} should be an object of class \code{"onearrow"}
  created by the command \code{\link{onearrow}}. 
}
\value{
  A window (class \code{"owin"}) enclosing the plotted graphics.
}
\examples{
  oa <- onearrow(cells[c(1, 42)])
  plot(oa)
  plot(oa, zap=TRUE, do.points=TRUE, col.head="pink", col="red")
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\seealso{
  \code{\link{onearrow}},
  \code{\link{yardstick}}
}
\keyword{spatial}
\keyword{hplot}
