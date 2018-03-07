\name{plot.tess}
\alias{plot.tess}
\title{Plot a Tessellation}
\description{
  Plots a tessellation, with optional labels for the tiles, and optional
  filled colour in each tile.
}
\usage{
  \method{plot}{tess}(x, \dots, main, add=FALSE,
                      show.all=!add,
                      border=NULL,
                      do.plot=TRUE,
                      do.labels=FALSE,
                      labels=tilenames(x), labelargs=list(),
                      do.col=FALSE,
                      values=marks(x), col=NULL, ribargs=list())
}
\arguments{
  \item{x}{Tessellation (object of class \code{"tess"}) to be plotted.}
  \item{\dots}{Arguments controlling the appearance of the plot.}
  \item{main}{Heading for the plot. A character string.}
  \item{add}{Logical. Determines whether the tessellation plot is added
    to the existing plot.
  }
  \item{show.all}{
    Logical value indicating whether to plot everything
    including the main title and the observation window of \code{x}.
  }
  \item{border}{
    Colour of the tile boundaries. A character string or other value
    specifying a single colour. 
    Ignored for pixel tessellations.
  }
  \item{do.plot}{
    Logical value indicating whether to actually perform the plot.
  }
  \item{do.labels}{
    Logical value indicating whether to show a text label for each tile
    of the tessellation.
  }
  \item{labels}{Character vector of labels for the tiles.}
  \item{labelargs}{
    List of arguments passed to
    \code{\link{text.default}} to control display of the text labels.
  }
  \item{do.col}{
    Logical value indicating whether tiles should be filled with
    colour. Always \code{TRUE} for pixel tessellations.
  }
  \item{values}{
    A vector of numerical values (or a factor, or vector of character
    strings) that will be associated with each tile of the tessellation
    and which determine the colour of the tile. The default is
    the marks of \code{x}. If the tessellation is not marked,
    or if the argument \code{values=NULL} is given, the default is
    a factor giving the tile identifier.
  }
  \item{col}{
    A vector of colours for each of the \code{values},
    or a \code{\link{colourmap}} that maps these values to colours.
  }
  \item{ribargs}{
    List of additional arguments to control the plot of the colour map,
    if \code{do.col=TRUE}. See explanation in \code{\link{plot.im}}.
  }
}
\details{
  This is a method for the generic \code{\link{plot}} function
  for the class \code{"tess"} of tessellations (see \code{\link{tess}}).

  The arguments \code{\dots} control the appearance of the plot.
  They are passed to \code{\link{segments}},
  \code{\link{plot.owin}} or \code{\link{plot.im}}, depending on the
  type of tessellation. 
}
\value{
  (Invisible) window of class \code{"owin"} specifying a bounding box
  for the plot, or an object of class \code{"colourmap"} specifying the
  colour map. (In the latter case, the bounding box information is
  available as an attribute, and can be extracted
  using \code{\link{as.owin}}.)
}
\seealso{
  \code{\link{tess}}
}
\examples{
  A <- tess(xgrid=0:4,ygrid=0:4)
  plot(A, border="blue", lwd=2, lty=2)
  plot(A, do.col=TRUE, values=runif(16), border="white")
  B <- A[c(1, 2, 5, 7, 9)]
  plot(B, hatch=TRUE)
  v <- as.im(function(x,y){factor(round(5 * (x^2 + y^2)))}, W=owin())
  levels(v) <- letters[seq(length(levels(v)))]
  E <- tess(image=v)
  plot(E)
  plot(E, col=rainbow(11), ribargs=list(las=1))
}

\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{hplot}
