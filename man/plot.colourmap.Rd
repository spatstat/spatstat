\name{plot.colourmap}
\alias{plot.colourmap}
\title{Plot a Colour Map}
\description{
  Displays a colour map as a colour ribbon
}
\usage{
\method{plot}{colourmap}(x, ...,
       main, xlim = NULL, ylim = NULL, vertical = FALSE, axis = TRUE,
       labelmap=NULL)
}
\arguments{
  \item{x}{Colour map to be plotted. An object of class \code{"colourmap"}.}
  \item{\dots}{
    Graphical arguments passed to \code{\link{image.default}}
    or \code{\link{axis}}.
  }
  \item{main}{Main title for plot. A character string.}
  \item{xlim}{
    Optional range of \code{x} values for the location of the
    colour ribbon.
  }
  \item{ylim}{
    Optional range of \code{y} values for the location of the
    colour ribbon.
  }
  \item{vertical}{Logical flag determining whether the colour ribbon
    is plotted as a horizontal strip (\code{FALSE}) or a vertical strip
    (\code{TRUE}).}
  \item{axis}{Logical flag determining whether an axis should be plotted
    showing the numerical values that are mapped to the colours.
  }
  \item{labelmap}{
    Function. If this is present, then the labels on the plot,
    which indicate the input values corresponding to particular colours,
    will be transformed by \code{labelmap} before being displayed
    on the plot. Typically used to simplify or shorten the
    labels on the plot.
  }
}
\details{
  This is the plot method for the class \code{"colourmap"}.
  An object of this class
  (created by the function \code{\link{colourmap}})
  represents a colour map or
  colour lookup table associating colours with each data value.

  The command \code{plot.colourmap} displays the colour map as a colour
  ribbon. This plot can be useful on its own to inspect the colour map.

  To annotate an existing plot with an explanatory colour ribbon,
  specify \code{add=TRUE} and use the arguments \code{xlim}
  and/or \code{ylim} to control the physical position of the ribbon
  on the plot.
}
\value{
  None.
}
\seealso{\code{\link{colourmap}}}
\examples{
  co <- colourmap(rainbow(100), breaks=seq(-1,1,length=101))
  plot(co)
  plot(co, vertical=TRUE)
  ca <- colourmap(rainbow(8), inputs=letters[1:8])
  plot(ca, col.ticks="pink")
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{color}
\keyword{hplot}
