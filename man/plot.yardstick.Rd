\name{plot.yardstick}
\alias{plot.yardstick}
\title{
  Plot a Yardstick or Scale Bar
}
\description{
  Plots an object of class \code{"yardstick"}.
}
\usage{
 \method{plot}{yardstick}(x, \dots,
        angle = 20, frac = 1/8, 
        split = FALSE, shrink = 1/4,
        pos = NULL,
        txt.args=list(),
        txt.shift=c(0,0),
        do.plot = TRUE)
}
\arguments{
  \item{x}{
    Object of class \code{"yardstick"} to be plotted.
    This object is created by the command \code{\link{yardstick}}.
  }
  \item{\dots}{
    Additional graphics arguments passed to
    \code{\link[graphics]{segments}} to control the appearance of the line.
  }
  \item{angle}{
    Angle between the arrows and the line segment, in degrees. 
  }
  \item{frac}{
    Length of arrow as a fraction of total length of the line segment.
  }
  \item{split}{
    Logical. If \code{TRUE}, then the line will be broken in the
    middle, and the text will be placed in this gap. If \code{FALSE},
    the line will be unbroken, and the text will be placed beside
    the line.
  }
  \item{shrink}{
    Fraction of total length to be removed from the middle of the
    line segment, if \code{split=TRUE}.
  }
  \item{pos}{
    Integer (passed to \code{\link[graphics]{text}}) determining the
    position of the annotation text relative to the line segment,
    if \code{split=FALSE}. Values of 1, 2, 3 and 4
    indicate positions below, to the left of,
    above and to the right of the line, respectively.
  }
  \item{txt.args}{
    Optional list of additional arguments passed to
    \code{\link[graphics]{text}} controlling the appearance of the text.
    Examples include \code{adj}, \code{srt}, \code{col}, \code{cex},
    \code{font}.
  }
  \item{txt.shift}{
    Optional numeric vector of length 2 specifying displacement
    of the text position relative to the centre of the yardstick.
  }
  \item{do.plot}{
    Logical. Whether to actually perform the plot (\code{do.plot=TRUE}).
  }
}
\details{
  A yardstick or scale bar is a line segment, drawn on any spatial
  graphics display, indicating the scale of the plot. 

  The argument \code{x} should be an object of class \code{"yardstick"}
  created by the command \code{\link{yardstick}}. 
}
\value{
  A window (class \code{"owin"}) enclosing the plotted graphics.
}
\examples{
  plot(owin(), main="Yardsticks")
  ys <- yardstick(as.psp(list(xmid=0.5, ymid=0.1, length=0.4, angle=0),
                         window=owin(c(0.2, 0.8), c(0, 0.2))),
                  txt="1 km")
  plot(ys)
  ys <- shift(ys, c(0, 0.3))
  plot(ys, angle=90, frac=0.08)
  ys <- shift(ys, c(0, 0.3))
  plot(ys, split=TRUE)
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\seealso{
  \code{\link{yardstick}}
}
\keyword{spatial}
\keyword{hplot}

