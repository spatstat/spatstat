\name{plot.lintess}
\alias{plot.lintess}
\title{
  Plot a Tessellation on a Linear Network
}
\description{
  Plot a tessellation or division of a linear network into tiles.
}
\usage{
\method{plot}{lintess}(x, \dots,
        main, add = FALSE, style = c("segments", "image"), col = NULL)
}
\arguments{
  \item{x}{
    Tessellation on a linear network
    (object of class \code{"lintess"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[graphics]{segments}} (if
    \code{style="segments"})
    or to \code{\link{plot.im}} (if \code{style="image"})
    to control the plot.
  }
  \item{main}{
    Optional main title for the plot.
  }
  \item{add}{
    Logical value indicating whether the plot is to be added to an
    existing plot.
  }
  \item{style}{
    Character string (partially matched) indicating whether to
    plot the tiles of the tessellation using \code{\link[graphics]{segments}}
    or to convert the tessellation to a pixel image and use
    \code{\link{plot.im}}.
  }
  \item{col}{
    Vector of colours, or colour map, determining the colours used
    to plot the different tiles of the tessellation.
  }
}
\details{
  A tessellation on a linear network \code{L} is a partition of the
  network into non-overlapping pieces (tiles). Each tile consists of one
  or more line segments which are subsets of the line segments making up
  the network. A tile can consist of several disjoint pieces.

  This function plots the tessellation on the current device.
  It is a method for the generic \code{plot}.
  
  If \code{style="segments"}, each tile is plotted using
  \code{\link[graphics]{segments}}. Colours distinguish the different
  tiles. 

  If \code{style="image"}, the tessellation is converted to a pixel
  image, and plotted using \code{\link{plot.im}}.
}
\value{
  (Invisible) colour map.
}
\author{
  \adrian
}
\seealso{
  \code{\link{lintess}}
}
\examples{
   X <- runiflpp(7, simplenet)
   Z <- divide.linnet(X)
   plot(Z, main="tessellation on network")
   points(as.ppp(X))
}
\keyword{spatial}
\keyword{hplot}
