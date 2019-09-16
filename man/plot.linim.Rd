\name{plot.linim}
\alias{plot.linim}
\title{
  Plot Pixel Image on Linear Network
}
\description{
  Given a pixel image on a linear network,
  the pixel values are displayed
  either as colours or as line widths.
}
\usage{
\method{plot}{linim}(x, ..., style = c("colour", "width"),
             scale, adjust = 1,
             negative.args = list(col=2),
             legend=TRUE,
             leg.side=c("right", "left", "bottom", "top"),
             leg.sep=0.1,
             leg.wid=0.1,
             leg.args=list(),
             leg.scale=1,
             zlim,
             box=FALSE,
             do.plot=TRUE)
}
\arguments{
  \item{x}{
    The pixel image to be plotted. An object of class \code{"linim"}.
  }
  \item{\dots}{
    Extra graphical parameters, passed to \code{\link{plot.im}}
    if \code{style="colour"}, or to \code{\link[graphics]{polygon}}
    if \code{style="width"}.
  }
  \item{style}{
    Character string (partially matched)
    specifying the type of plot. See Details.
  }
  \item{scale}{
    Physical scale factor for representing the pixel values as
    line widths. 
  }
  \item{adjust}{
    Adjustment factor for the default scale.
  }
  \item{negative.args}{
    A list of arguments to be passed to \code{\link[graphics]{polygon}}
    specifying how to plot negative values of \code{x}
    when \code{style="width"}.
  }
  \item{legend}{
    Logical value indicating whether to plot a legend
    (colour ribbon or scale bar).
  }
  \item{leg.side}{
    Character string (partially matched) indicating where to display the legend
    relative to the main image.
  }
  \item{leg.sep}{
    Factor controlling the space between the legend and the image.
  }
  \item{leg.wid}{
    Factor controlling the width of the legend.
  }
  \item{leg.scale}{
    Rescaling factor for annotations on the legend.
    The values on the numerical scale
    printed beside the legend will be multiplied by this rescaling factor.
  }
  \item{leg.args}{
    List of additional arguments passed to
    \code{\link[graphics]{image.default}},
    \code{\link[graphics]{axis}}
    or \code{\link[graphics]{text.default}} 
    to control the display of the legend.
    These may override the \code{\dots} arguments. 
  }
  \item{zlim}{
    The range of numerical values that should be mapped.
    A numeric vector of length 2.
    Defaults to the range of values of \code{x}.
  }
  \item{box}{
    Logical value indicating whether to draw a bounding box.
  }
  \item{do.plot}{
    Logical value indicating whether to actually perform the plot.
  }
}
\details{
  This is the \code{plot} method for objects 
  of class \code{"linim"}. Such an object represents
  a pixel image defined on a linear network.

  If \code{style="colour"} (the default) then
  the pixel values of \code{x} are plotted as colours,
  using \code{\link{plot.im}}.
  The mapping from pixel values to colours is determined by any
  additional arguments \code{\dots} which are passed to \code{\link{plot.im}}.

  If \code{style="width"} then
  the pixel values of \code{x} are used to determine the widths of
  thick lines centred on the line segments of the linear network.
  The mapping from pixel values to line widths is determined by
  the arguments \code{scale} and \code{adjust}.
  The plotting of colours and borders of the lines is controlled by
  the additional arguments \code{\dots} which are passed to
  \code{\link[graphics]{polygon}}. A different set of colours and borders can be
  assigned to negative pixel values by passing a list of arguments in
  \code{negative.args} as shown in the Examples.

  A legend is displayed alongside the plot
  if \code{legend=TRUE} (the default). The legend
  displays the relationship between pixel values and colours
  (if \code{style="colour"}) or between pixel values and
  line widths (if \code{style="width"}).

  The plotting of the legend itself is controlled by the arguments
  \code{leg.side}, \code{leg.sep}, \code{leg.wid}, \code{leg.scale}
  and the list of arguments \code{leg.args}, which are described above.
  If \code{style="colour"}, these arguments are mapped to
  the arguments \code{ribside}, \code{ribsep}, \code{ribwid},
  \code{ribscale} and \code{ribargs} respectively, which are passed
  to \code{\link{plot.im}}.
}
\value{
  If \code{style="colour"}, the result is
  an object of class \code{"colourmap"} specifying the colour map used.
  If \code{style="width"}, the result is 
  a numeric value \code{v} giving the physical scale:
  one unit of pixel value is represented as \code{v} physical units on the plot.

  The result also has an attribute \code{"bbox"} giving a bounding box
  for the plot. The bounding box includes the ribbon or scale bar, if present,
  but not the main title.
}
\author{
  \adrian
}
\seealso{
  \code{\link{linim}},
  \code{\link{plot.im}},
  \code{\link[graphics]{polygon}}
}
\references{
  Ang, Q.W., Baddeley, A. and Nair, G. (2012)
  Geometrically corrected second-order analysis of 
  events on a linear network, with applications to
  ecology and criminology.
  \emph{Scandinavian Journal of Statistics} \bold{39}, 591--617.
}
\examples{
  X <- linfun(function(x,y,seg,tp){y^2+x}, simplenet)
  X <- as.linim(X)
  
  plot(X)
  plot(X, style="width", main="Width proportional to function value")

  # signed values
  f <- linfun(function(x,y,seg,tp){y-x}, simplenet)
  plot(f, style="w", main="Negative values in red")

  plot(f, style="w", negative.args=list(density=10),
       main="Negative values are hatched")
}
\keyword{spatial}
