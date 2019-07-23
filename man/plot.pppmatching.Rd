\name{plot.pppmatching}
\alias{plot.pppmatching}
\title{
  Plot a Point Matching
}
\description{
  Plot an object of class \code{"pppmatching"} which represents
  a matching of two planar point patterns.
}
\usage{
\method{plot}{pppmatching}(x, addmatch = NULL, main = NULL, \dots, adjust = 1)
}
\arguments{
  \item{x}{
    Point pattern matching object (class \code{"pppmatching"})
    to be plotted.
  }
  \item{addmatch}{
    Optional. A matrix indicating additional pairs of points that
    should be matched. See Details.
  }
  \item{main}{
    Main title for the plot. 
  }
  \item{\dots}{
    Additional arguments passed to other plot methods.
  }
  \item{adjust}{
    Adjustment factor for the widths of line segments.
    A positive number.
  }
}
\details{
  The object \code{x} represents a matching found between
  two point patterns \code{X} and \code{Y}. The matching may be
  incomplete. See \code{\link{pppmatching.object}} for further
  description.

  This function plots the matching by drawing the two point patterns
  \code{X} and \code{Y} as red and blue dots respectively, and drawing
  line segments between each pair of matched points. The width of the
  line segments is proportional to the strength of matching. The
  proportionality constant can be adjusted using the argument
  \code{adjust}.

  Additional graphics arguments \code{\dots} control the
  plotting of the window (and are passed to \code{\link{plot.owin}})
  and the plotting of the line segments 
  (and are passed to \code{\link{plot.psp}}, \code{\link{plot.linim}}
  and ultimately to the base graphics function
  \code{\link[graphics]{polygon}}).

  The argument \code{addmatch} is for use mainly by developers to
  study algorithms which update the matching.
  If \code{addmatch} is given, it should be a matrix with dimensions
  \code{npoints(X) * npoints(Y)}. If \code{addmatch[i,j] > 0} then
  a light grey line segment will be drawn between \code{X[i]} and \code{Y[j}.
}
\value{
  Null.
}
\author{
  Dominic Schuhmacher and \adrian.
}
\seealso{
  \code{\link{pppmatching.object}}
}
\examples{
  X <- runifpoint(7)
  Y <- runifpoint(7)
  am <- r2dtable(1, rep(10,7), rep(10,7))[[1]]/10
  m2 <- pppmatching(X, Y, am)
  plot(m2, adjust=0.3)
}
\keyword{spatial}
\keyword{hplot}
