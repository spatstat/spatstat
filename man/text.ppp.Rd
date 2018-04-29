\name{text.ppp}
\alias{text.ppp}
\alias{text.lpp}
\alias{text.psp}
\title{
  Add Text Labels to Spatial Pattern 
}
\description{
  Plots a text label at the location of each point in a spatial point
  pattern, or each object in a spatial pattern of objects.
}
\usage{
\method{text}{ppp}(x, \dots)

\method{text}{lpp}(x, \dots)

\method{text}{psp}(x, \dots)
}
\arguments{
  \item{x}{
    A spatial point pattern (object of class \code{"ppp"}),
    a point pattern on a linear network (class \code{"lpp"})
    or a spatial pattern of line segments (class \code{"psp"}).
  }
  \item{\dots}{
   Additional arguments passed to \code{\link[graphics]{text.default}}.
  }
}
\details{
  These functions are methods for the generic \code{\link{text}}.
  A text label is added to the existing plot,
  at the location of each point in the point pattern \code{x},
  or near the location of the midpoint of each segment
  in the segment pattern \code{x}.

  Additional arguments \code{\dots} are passed to
  \code{\link[graphics]{text.default}} and may be used to control the
  placement of the labels relative to the point locations, and
  the size and colour of the labels.

  By default, the labels are the serial numbers 1 to \code{n}, where
  \code{n} is the number of points or segments in \code{x}. This can
  be changed by specifying the argument \code{labels}, which should be
  a vector of length \code{n}.
}
\value{
  Null.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[graphics]{text.default}}
}
\examples{
  plot(cells)
  text(cells, pos=2)

  plot(Frame(cells))
  text(cells, cex=1.5)

  S <- as.psp(simplenet)
  plot(S)
  text(S)

  X <- runiflpp(5, simplenet)
  plot(X)
  text(X, pos=2, col="blue")
}
\keyword{spatial}
\keyword{hplot}
