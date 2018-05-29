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
                      values=marks(x),
                      multiplot=TRUE, 
                      col=NULL, ribargs=list())
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
  \item{multiplot}{
    Logical value giving permission to display multiple plot panels.
    This applies when \code{do.col=TRUE} and \code{ncol(values) > 1}.
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

  The window of the tessellation is plotted, and then the tiles of the
  tessellation are plotted in their correct positions in the window.

  Rectangular or polygonal tiles are plotted individually
  using \code{\link{plot.owin}}, while a
  tessellation represented by a pixel image
  is plotted using \code{\link{plot.im}}.
  The arguments \code{\dots} control the appearance of the plot,
  and are passed to \code{\link{segments}},
  \code{\link{plot.owin}} or \code{\link{plot.im}} as appropriate.

  If \code{do.col=TRUE}, then the tiles of the tessellation are filled
  with colours determined by the argument \code{values}.
  By default, these values are the marks associated with each of the tiles.
  If there is more than one column of marks or values,
  then the default behaviour (if \code{multiplot=TRUE})
  is to display several plot panels, one for
  each column of mark values. Then the arguments \code{\dots} are
  passed to \code{\link{plot.solist}} to determine the arrangement of
  the panels. 
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
  Rect <- tess(xgrid=0:4,ygrid=0:4)
  Diri <- dirichlet(runifpoint(7))
  plot(Diri)
  plot(Rect, border="blue", lwd=2, lty=2)
  plot(Rect, do.col=TRUE, border="white")
  plot(Rect, do.col=TRUE, values=runif(16), border="white")
  B <- Rect[c(1, 2, 5, 7, 9)]
  plot(B, hatch=TRUE)
  plot(Diri, do.col=TRUE)
  plot(Diri, do.col=TRUE, do.labels=TRUE, labelargs=list(col="white"),
             ribbon=FALSE)
  v <- as.im(function(x,y){factor(round(5 * (x^2 + y^2)))}, W=owin())
  levels(v) <- letters[seq(length(levels(v)))]
  Img <- tess(image=v)
  plot(Img)
  plot(Img, col=rainbow(11), ribargs=list(las=1))
  a <- tile.areas(Diri)
  marks(Diri) <- data.frame(area=a, random=runif(7, max=max(a)))
  plot(Diri, do.col=TRUE, equal.ribbon=TRUE)
}

\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{hplot}
