\name{textureplot}
\alias{textureplot}
\title{
  Plot Image or Tessellation Using Texture Fill
}
\description{
  For a factor-valued pixel image, this command plots each level
  of the factor using a different texture.
  For a tessellation, each tile is plotted using a different texture.
}
\usage{
textureplot(x, \dots,
            main, add=FALSE, clipwin=NULL, do.plot = TRUE,
            border=NULL, col = NULL, lwd = NULL, lty = NULL, spacing = NULL,
            textures=1:8,
            legend=TRUE,
            leg.side=c("right", "left", "bottom", "top"),
            legsep=0.1, legwid=0.2)
}
\arguments{
  \item{x}{
    A tessellation (object of class \code{"tess"} or something
    acceptable to \code{\link{as.tess}}) with at most 8 tiles,
    or a pixel image (object of class \code{"im"} or something
    acceptable to \code{\link{as.im}}) whose pixel values
    are a \code{factor} with at most 8 levels.
  }
  \item{\dots}{
    Other arguments passed to \code{\link{add.texture}}.
  }
  \item{main}{
    Character string giving a main title for the plot.
  }
  \item{add}{
    Logical value indicating whether to draw on
    the current plot (\code{add=TRUE}) or to initialise a new plot
    (\code{add=FALSE}).
  }
  \item{clipwin}{
    Optional. A window (object of class \code{"owin"}).
    Only this subset of the image will be displayed.
  }
  \item{do.plot}{
    Logical. Whether to actually do the plot.
  }
  \item{border}{
    Colour for drawing the boundaries between the different regions.
    The default (\code{border=NULL}) means to use \code{par("fg")}.
    Use \code{border=NA} to omit borders.
  }
  \item{col}{
    Numeric value or vector giving the colour or colours in which
    the textures should be plotted.
  }
  \item{lwd}{
    Numeric value or vector giving the line width or widths to be used.
  }
  \item{lty}{
    Numeric value or vector giving the line type or types to be used.
  }
  \item{spacing}{
    Numeric value or vector giving the \code{spacing} parameter
    for the textures.
  }
  \item{textures}{
    Textures to be used for each level. Either a texture map (object of
    class \code{"texturemap"}) or a vector of integer codes
    (to be interpreted by \code{\link{add.texture}}).
  }
  \item{legend}{
    Logical. Whether to display an explanatory legend.
  }
  \item{leg.side}{Position of legend relative to main plot.}
  \item{legsep}{
    Separation between legend and main plot, as a fraction
    of the shortest side length of the main plot.
  }
  \item{legwid}{
    Width (if vertical) or height (if horizontal) of the legend
    as a fraction of the shortest side length of the main plot.
  }
}
\details{
  If \code{x} is a tessellation, then each tile of the tessellation is
  plotted and filled with a texture using \link{add.texture}.

  If \code{x} is a factor-valued pixel image, then
  for each level of the factor, the algorithm finds the region where the image
  takes this value, and fills the region with a texture using
  \code{\link{add.texture}}.
}
\value{
  (Invisible) A texture map (object of class \code{"texturemap"})
  associating a texture with each level of the factor.
}
\author{
Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}


Rolf Turner \email{r.turner@auckland.ac.nz}

and Ege Rubak \email{rubak@math.aau.dk}

}
\seealso{
  \code{\link{im}}, 
  \code{\link{plot.im}},
  \code{\link{add.texture}}.
}
\examples{
  Z <- cut(setcov(owin()), 3, labels=c("Lo", "Med", "Hi"))
  textureplot(Z)
  textureplot(dirichlet(runifpoint(6)))
}
\keyword{spatial}
\keyword{hplot}
