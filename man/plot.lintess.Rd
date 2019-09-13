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
        main, add = FALSE,
        style = c("colour", "width", "image"),
        col = NULL, values=marks(x),
        ribbon=TRUE, ribargs=list(), multiplot=TRUE, do.plot=TRUE)
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
    Character string (partially matched) specifying the type of plot.
    If \code{style="colour"} (the default), tiles are plotted
    using \code{\link[graphics]{segments}} using colours to distinguish
    the different tiles or values. If \code{style="width"},
    tiles are plotted using \code{\link[graphics]{segments}}
    using different segment widths to distinguish
    the different tiles or values. If \code{style="image"}, the
    tessellation is converted to a pixel image and plotted by
    \code{\link{plot.im}}.
  }
  \item{col}{
    Vector of colours, or colour map, determining the colours used
    to plot the different tiles of the tessellation.
  }
  \item{values}{
    Values associated with each tile of the tessellation,
    used to determine the colours or widths.
    A vector with one entry for each tile, or a data frame
    with one row for each tile.
    The default is \code{marks(x)}, or if that is null, then
    \code{tilenames(x)}.
  }
  \item{ribbon}{
    Logical value specifying whether to print an explanatory legend
    for the colour map or width map.
  }
  \item{ribargs}{
    Arguments passed to \code{\link{plot.colourmap}} controlling the
    display of the colour map legend.
  }
  \item{multiplot}{
    Logical value determining what should happen if \code{marks(x)}
    has more than one column. If \code{multiplot=TRUE} (the default),
    several plot panels will be generated, one panel for each column
    of marks. If \code{multiplot=FALSE}, the first column of marks will
    be selected.
  }
  \item{do.plot}{
    Logical value specifying whether to actually generate the plot
    (\code{do.plot=TRUE}, the default) or just to compute the colour map
    and return it (\code{do.plot=FALSE}).
  }
}
\details{
  A tessellation on a linear network \code{L} is a partition of the
  network into non-overlapping pieces (tiles). Each tile consists of one
  or more line segments which are subsets of the line segments making up
  the network. A tile can consist of several disjoint pieces.

  This function plots the tessellation on the current device.
  It is a method for the generic \code{plot}.

  If \code{style="colour"}, each tile is plotted using
  \code{\link[graphics]{segments}}, drawing segments of different colours. 
  
  If \code{style="width"}, each tile is plotted using
  \code{\link[graphics]{segments}}, drawing segments of different widths.
  
  If \code{style="image"}, the tessellation is converted to a pixel
  image, and plotted as a colour image using \code{\link{plot.im}}.

  The colours or widths are determined by the \code{values}
  associated with each tile of the tessellation. If \code{values} is
  missing, the default is to use the marks of the tessellation, or if
  there are no marks, the names of the tiles.
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
   plot(Z, main="tessellation on network",
           values=1:nobjects(Z), style="w")
}
\keyword{spatial}
\keyword{hplot}
