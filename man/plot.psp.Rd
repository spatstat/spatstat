\name{plot.psp}
\alias{plot.psp}
\title{plot a Spatial Line Segment Pattern}
\description{
  Plot a two-dimensional line segment pattern
}
\usage{
  \method{plot}{psp}(x, \dots, main, add=FALSE,
                     show.all=!add, show.window=show.all, do.plot=TRUE,
                     which.marks=1,
                     style=c("colour", "width", "none"),
                     col=NULL, 
                     ribbon=show.all,
                     ribsep=0.15, ribwid=0.05, ribn=1024,
                     scale=NULL, adjust=1,
                     legend=TRUE,
                     leg.side=c("right", "left", "bottom", "top"),
                     leg.sep=0.1,
                     leg.wid=0.1,
                     leg.args=list(),
                     leg.scale=1,
                     negative.args=list(col=2))
}
\arguments{
  \item{x}{
    The line segment pattern to be plotted.
    An object of class \code{"psp"},
    or data which can be converted into 
    this format by \code{\link{as.psp}()}.
  }
  \item{\dots}{
    extra arguments that will be passed to the plotting functions
    \code{\link{segments}} (to plot the segments)
    and 
    \code{\link{plot.owin}} (to plot the observation window).
  }
  \item{main}{
    Character string giving a title for the plot.
  }
  \item{add}{
    Logical. If \code{TRUE}, the current plot is not erased;
    the segments are plotted on top of the
    current plot, and the window is not plotted (by default).
  }
  \item{show.all}{
    Logical value specifying whether to plot everything
    including the window, main title, and colour ribbon.
  }
  \item{show.window}{
    Logical value specifying whether to plot the window.
  }
  \item{do.plot}{
    Logical value indicating whether to actually perform the plot.
  }
  \item{which.marks}{
    Index determining which column of marks to use,
    if the marks of \code{x} are a data frame.
    A character string or an integer.
    Defaults to \code{1} indicating the first column of marks.
  }
  \item{style}{
    Character string specifying how to represent the
    mark value of each segment. If \code{style="colour"} (the default)
    segments are coloured according to their mark value.
    If \code{style="width"}, segments are drawn with a width
    proportional to their mark value.
    If \code{style="none"} the mark values are ignored.
  }
  \item{col}{
    Colour information.
    If \code{style="width"} or \code{style="none"},
    then \code{col} should be a single value,
    interpretable as a colour; the line segments will be plotted
    using this colour.
    If \code{style="colour"} and \code{x} has marks,
    then the mark values will be mapped
    to colours using the information in \code{col}, which should be a 
    colour map (object of class \code{"colourmap"})
    or a vector of colour values.
  }
  \item{ribbon}{
    Logical value indicating whether to display a ribbon
    showing the colour map (in which mark values are associated with
    colours) when \code{style="colour"}.
  }
  \item{ribsep}{
    Factor controlling the space between the colour ribbon and the image.
  }
  \item{ribwid}{
    Factor controlling the width of the colour ribbon.
  }
  \item{ribn}{
    Number of different values to display in the colour ribbon.
  }
  \item{scale}{
    Optional. Physical scale for representing the mark values of
    \code{x} as physical widths on the plot, when \code{style="width"}.
    There is a sensible default.
  }
  \item{adjust}{
    Optional adjustment factor for \code{scale}.
  }
  \item{legend}{
    Logical value indicating whether to display a legend
    showing the width map (in which mark values are associated with
    segment widths) when \code{style="width"}.
  }
  \item{leg.side}{
    Character string (partially matched) specifying where the legend
    should be plotted,
    when \code{style="width"}.
  }
  \item{leg.sep}{
    Factor controlling the space between the legend and the main plot,
    when \code{style="width"}.
  }
  \item{leg.wid}{
    Factor controlling the width of the legend,
    when \code{style="width"}.
  }
  \item{leg.args}{
    Optional list of additional arguments passed to
    \code{\link[graphics]{axis}} and
    \code{\link[graphics]{text.default}}
    controlling the appearance of the legend,
    when \code{style="width"}.
  }
  \item{leg.scale}{
    Rescaling factor for labels, when \code{style="width"}.
    The values on the numerical scale
    printed beside the legend will be multiplied by this rescaling factor.
  }
  \item{negative.args}{
    Optional list of arguments to \code{\link[graphics]{polygon}}
    to be used when the mark values are negative. 
  }
}
\value{
  If \code{style="colour"}, the result is a \code{\link{colourmap}}
  object specifying the association between
  marks and colours, if any.
  
  If \code{style="width"}, the result is a numeric value giving the
  scaling between the mark values and the physical widths.

  In all cases, the return value also has an attribute
  \code{"bbox"} giving a bounding box for the plot.
}
\details{
  This is the \code{plot} method for 
  line segment pattern datasets (of class \code{"psp"},
  see \code{\link{psp.object}}).
  It plots both the observation window \code{Window(x)}
  and the line segments themselves.
  
  Plotting of the window \code{Window(x)} is performed by
  \code{\link{plot.owin}}. This plot may be modified 
  through the \code{...} arguments. 

  Plotting of the segments themselves is performed
  by the standard R function \code{\link{segments}}.
  Its plotting behaviour may also be modified through the \code{...}
  arguments.

  There are three different styles of plotting
  which apply when the segments have marks
  (i.e. when \code{marks(x)} is not null):
  \describe{
    \item{\code{style="colour"} (the default):}{
      Segments are plotted with different colours depending on their
      mark values.
      The colour map, associating mark values with colours,
      is determined by the argument \code{col}. The colour map will be
      displayed as a vertical colour ribbon to the right of the plot, if
      \code{ribbon=TRUE} (the default).
    }
    \item{\code{style="width"}:}{
      Segments are plotted with different widths depending on their
      mark values. 
      The expanded segments are plotted using the base graphics
      function \code{\link[graphics]{polygon}}.
      The width map, associating mark values with line widths,
      can be specified by giving the physical scale factor \code{scale}.
      There is a sensible default scale, which can be adjusted using
      the adjustment factor \code{adjust}.
      The width map will be displayed as a vertical stack of lines
      to the right of the plot, if \code{legend=TRUE} (the default).
    }
    \item{\code{style="none"}:}{
      Mark information is ignored and the
      segments are plotted as thin lines using \code{\link[graphics]{segments}}.
    }
  }
  If \code{marks(x)} is a data frame, the default is to use the first
  column of \code{marks(x)} to determine the colours or widths.
  To specify another
  column, use the argument \code{which.marks}.
}
\seealso{
  \code{\link{psp.object}},
  \code{\link{plot}},
  \code{\link{par}},
  \code{\link{plot.owin}},
  \code{\link{text.psp}},
  \code{\link{symbols}}
}
\examples{
  X <- psp(runif(20), runif(20), runif(20), runif(20), window=owin())
  plot(X)
  plot(X, lwd=3)
  lettuce <- sample(letters[1:4], 20, replace=TRUE)
  marks(X) <- data.frame(A=1:20, B=factor(lettuce))
  plot(X)
  plot(X, which.marks="B")
  plot(X, style="width", col="grey")
}

\author{
  \adrian
  and \rolf
}
\keyword{spatial}
 
 \keyword{hplot}
