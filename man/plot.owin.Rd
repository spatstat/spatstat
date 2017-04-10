\name{plot.owin}
\alias{plot.owin}
\title{Plot a Spatial Window}
\description{
  Plot a two-dimensional window of observation for a spatial point pattern
}
\usage{
 \method{plot}{owin}(x, main, add=FALSE, \dots, box, edge=0.04,
                      type=c("w","n"), show.all=!add, 
                      hatch=FALSE,
                      hatchargs=list(), 
                      invert=FALSE, do.plot=TRUE,
                      claim.title.space=FALSE)
}
\arguments{
  \item{x}{
    The window to be plotted.
    An object of class \code{\link{owin}},
    or data which can be converted into 
    this format by \code{\link{as.owin}()}.
  }
  \item{main}{
    text to be displayed as a title above the plot.
  }
  \item{add}{
    logical flag: if \code{TRUE}, draw the window in 
    the current plot; if \code{FALSE}, generate a new plot.
  }
  \item{\dots}{
    extra arguments controlling the appearance of the plot.
    These arguments are passed to \code{\link[graphics]{polygon}}
    if \code{x} is a
    polygonal or rectangular window, or passed to
    \code{\link[graphics]{image.default}}
    if \code{x} is a binary mask. See Details.
  }
  \item{box}{
    logical flag; if \code{TRUE}, plot the enclosing rectangular box
  }
  \item{edge}{
    nonnegative number; the plotting region will have coordinate limits
    that are \code{1 + edge}
    times as large as the limits of the rectangular box
    that encloses the pattern.
  }
  \item{type}{
    Type of plot: either \code{"w"} or \code{"n"}.
    If \code{type="w"} (the default), the window is plotted.
    If \code{type="n"} and \code{add=TRUE}, a new plot is initialised
    and the coordinate system is established, but nothing is drawn.
  }
  \item{show.all}{
    Logical value indicating whether to plot everything
    including the main title.
  }
  \item{hatch}{
    logical flag; if \code{TRUE}, the interior of the window will
    be shaded by texture, such as a grid of parallel lines.
  }
  \item{hatchargs}{
    List of arguments passed to \code{\link{add.texture}}
    to control the texture shading when \code{hatch=TRUE}. 
  }
  \item{invert}{
    logical flag; when the window is a binary pixel mask,
    the mask colours will be inverted if \code{invert=TRUE}.
  }
  \item{do.plot}{
    Logical value indicating whether to actually perform the plot.
  }
  \item{claim.title.space}{
    Logical value indicating whether extra space for the main title
    should be allocated when declaring the plot dimensions.
    Should be set to \code{FALSE} under normal conditions.
  }
}
\value{
  none.
}
\details{
  This is the \code{plot} method for the class \code{\link{owin}}.
  The action is to plot the boundary of the window on the current plot device,
  using equal scales on the \code{x} and \code{y} axes.

  If the window \code{x} is of type \code{"rectangle"} or \code{"polygonal"},
  the boundary of the window is plotted as a polygon or series of polygons.
  If \code{x} is of type \code{"mask"} 
  the discrete raster approximation of the window is displayed
  as a binary image (white inside the window, black outside).

  Graphical parameters controlling the display (e.g. setting the
  colours) may be passed directly via the \code{...} arguments,
  or indirectly reset using 
  \code{\link{spatstat.options}}.

  When \code{x} is of type \code{"rectangle"} or \code{"polygonal"}, it
  is plotted by the \R function \code{\link[graphics]{polygon}}. To control the
  appearance (colour, fill density, line density etc) of the polygon
  plot, determine the required argument of \code{\link[graphics]{polygon}} and
  pass it through \code{...} For example, to paint the interior of the
  polygon in red, use the argument \code{col="red"}. To draw the polygon
  edges in green, use \code{border="green"}. To suppress the drawing of
  polygon edges, use \code{border=NA}.

  When \code{x} is of type \code{"mask"}, it is plotted by
  \code{\link[graphics]{image.default}}. The appearance of the image plot
  can be controlled by passing arguments to
  \code{\link[graphics]{image.default}}
  through \code{...}. The default appearance can also be changed
  by setting the parameter \code{par.binary} of
  \code{\link{spatstat.options}}. 

  To zoom in (to view only a subset of the window at higher
  magnification), use the graphical arguments
  \code{xlim} and \code{ylim} to specify the desired rectangular field of
  view. (The actual field of view may be larger, depending on the
  graphics device).
}
\section{Notes on Filled Polygons with Holes}{
  The function \code{\link[graphics]{polygon}} can only handle
  polygons without holes. To plot polygons with holes in a solid colour,
  we have implemented two workarounds.
  
  \describe{
    \item{polypath function:}{
      The first workaround uses
      the relatively new function \code{\link[graphics]{polypath}} which
      \emph{does} have the capability to handle polygons with holes.
      However, not all graphics devices support
      \code{\link[graphics]{polypath}}.
      The older devices \code{\link{xfig}} and \code{\link{pictex}}
      do not support \code{\link[graphics]{polypath}}.
      On a Windows system, the default graphics device
#ifdef windows      
      \code{\link{windows}} 
#endif
#ifndef windows      
      \code{windows}
#endif
      supports \code{\link[graphics]{polypath}}.
#ifdef unix
      On a Linux system, the default graphics device
      \code{X11(type="Xlib")} does \emph{not} support
      \code{\link[graphics]{polypath}}
      but \code{X11(type="cairo")} does support it.
      See \code{\link{X11}} and the section on Cairo below.
#endif
    }
    \item{polygon decomposition:}{
      The other workaround involves decomposing the polygonal window
      into pieces which do not have holes. This code is experimental
      but works in all our test cases. If this code fails, a warning
      will be issued, and the filled colours will not be plotted.
    }
  }
}
#ifdef unix
\section{Cairo graphics on a Linux system}{
  Linux systems support
  the graphics device \code{X11(type="cairo")} (see \code{\link{X11}})
  provided the external library \pkg{cairo} is installed
  on the computer. See \code{www.cairographics.org}
  for instructions on obtaining and installing \pkg{cairo}.  After having
  installed \pkg{cairo} one needs to re-install \R from source so
  that it has \pkg{cairo} capabilites.  To check whether your
  current installation of R has \pkg{cairo} capabilities, type
  (in \R) \code{capabilities()["cairo"]}.
  The default type for \code{\link{X11}} is controlled by
  \code{\link[grDevices]{X11.options}}.
  You may find it convenient to
  make \pkg{cairo} the default, e.g. via your \code{.Rprofile}.
  The magic incantation to put into \code{.Rprofile} is
  \preformatted{
    setHook(packageEvent("graphics", "onLoad"),
    function(...) grDevices::X11.options(type="cairo"))
  }
}
#endif
\seealso{
  \code{\link{owin.object}},
  \code{\link{plot.ppp}},
  \code{\link[graphics]{polygon}},
  \code{\link[graphics]{image.default}},
  \code{\link{spatstat.options}}
}
\examples{
  # rectangular window
   plot(Window(nztrees))
   abline(v=148, lty=2)

  # polygonal window
  w <- Window(demopat)
  plot(w)
  plot(w, col="red", border="green", lwd=2)
  plot(w, hatch=TRUE, lwd=2)

  # binary mask
  we <- as.mask(w)
  plot(we)
  op <- spatstat.options(par.binary=list(col=grey(c(0.5,1))))
  plot(we)
  spatstat.options(op)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{hplot}
