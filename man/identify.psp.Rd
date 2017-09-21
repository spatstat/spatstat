\name{identify.psp}
\alias{identify.psp}
\title{Identify Segments in a Line Segment Pattern}
\description{
  If a line segment pattern is plotted in the graphics window,
  this function will find the segment which is nearest to
  the mouse position, and print its serial number.
}
\usage{
  \method{identify}{psp}(x, \dots, labels=seq_len(nsegments(x)), n=nsegments(x), plot=TRUE)
}
\arguments{
  \item{x}{
    A line segment pattern (object of class \code{"psp"}).
  }
  \item{labels}{
    Labels associated with the segments, to be plotted when the
    segments are identified. A character vector or numeric vector
    of length equal to the number of segments in \code{x}.
  }
  \item{n}{
    Maximum number of segments to be identified.
  }
  \item{plot}{
    Logical. Whether to plot the labels when a segment is identified.
  }
  \item{\dots}{
    Arguments passed to \code{\link[graphics]{text.default}}
    controlling the plotting of the labels.
  }
}
\value{
  Vector containing the serial numbers of the segments in the pattern
  \code{x} that were identified.
}
\details{
  This is a method for the generic function \code{\link[graphics]{identify}}
  for line segment pattern objects.

  The line segment pattern \code{x} should first be plotted
  using \code{\link{plot.psp}}. Then \code{identify(x)}
  reads the position of the graphics pointer each time the
  left mouse button is pressed.  It then finds 
  the segment in the pattern \code{x} that is closest to the mouse position.
  This segment's index will be returned as part of the value of the call.

  Each time a segment is identified,
  text will be displayed next to the point,
  showing its serial number (or the relevant entry of \code{labels}).
}
\seealso{
  \code{\link[graphics]{identify}},
  \code{\link{identify.ppp}}.
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{iplot}
