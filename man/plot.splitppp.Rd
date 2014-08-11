\name{plot.splitppp}
\alias{plot.splitppp}
\title{Plot a List of Point Patterns}
\description{
  Plots a list of point patterns.
}
\usage{
  \method{plot}{splitppp}(x, \dots, main, arrange=TRUE,
   nrows=NULL, ncols=NULL,
   main.panel=NULL, mar.panel=c(2,1,1,2),
   panel.begin=NULL, panel.end=NULL, panel.args=NULL, plotcommand="plot",
   adorn.left=NULL, adorn.right=NULL, adorn.top=NULL, adorn.bottom=NULL,
   adorn.size=0.2,
   equal.scales=FALSE)
}
\arguments{
  \item{x}{
    A named list of point patterns,
    typically obtained from \code{\link{split.ppp}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.ppp}} which control the
    appearance of each plot panel.
  }
  \item{main}{
    Overall heading for the plot.
  }
  \item{arrange}{
    Logical flag indicating whether to plot the point patterns
    side-by-side on a single page (\code{arrange=TRUE})
    or plot them individually in a succession of frames
    (\code{arrange=FALSE}).
  }
  \item{nrows,ncols}{
    Optional. The number of rows/columns in the plot layout
    (assuming \code{arrange=TRUE}).
    You can specify either or both of these numbers.
  }
  \item{main.panel}{
    Optional. A character string, or a vector of character strings,
    giving the headings for each of the point patterns. 
  }
  \item{mar.panel}{
    Optional value of the graphics parameter \code{mar} controlling the size
    of the margins outside each plot panel.
    See the help file for \code{\link{par}}.
  }
  \item{panel.begin,panel.end}{
    Optional. Functions
    that will be executed before and after each panel is
    plotted.
    See Details.
  }
  \item{panel.args}{
    Internal use only.
  }
  \item{plotcommand}{
    Optional.
    Character string containing the name of the command that should be
    executed to plot each panel. 
  }
  \item{adorn.left,adorn.right,adorn.top,adorn.bottom}{
    Optional. Functions (with no arguments) that will be executed to
    generate additional plots at the margins (left, right, top and/or
    bottom, respectively) of the array of plots.
  }
  \item{adorn.size}{
    Relative width (as a fraction of the other panels' widths)
    of the margin plots.
  }
  \item{equal.scales}{
    Logical value indicating whether the components
    should be plotted at (approximately) the same physical scale. 
  }
}
\value{
  Null.
}
\details{
  This is the \code{plot} method for the class \code{"splitppp"}.
  It is typically used to plot the result of the function
  \code{\link{split.ppp}} but it may also be used to plot any list
  of point patterns created by the user.

  The argument \code{x} should be a named list of point patterns
  (objects of class \code{"ppp"}, see \code{\link{ppp.object}}).
  Each of these point patterns will be plotted in turn
  using \code{\link{plot.ppp}}. 

  The arguments \code{panel.begin} and \code{panel.end}
  may be functions that will be executed before and after each panel is plotted.
  They will be called as \code{panel.begin(i, y, main=main.panel[i])}
  and \code{panel.end(i, y, add=TRUE)}.

  Alternatively, \code{panel.begin} and \code{panel.end} may be objects
  of some class that can be plotted
  with the generic \code{plot} command. They will be plotted before and
  after each panel is plotted.

  If \code{equal.scales=FALSE} (the default), then the 
  plot panels will have equal height on the plot device
  (unless there is only one column of panels, in which case
  they will have equal width on the plot device). This means that the
  objects are plotted at different physical scales, by default.

  If \code{equal.scales=TRUE}, then the dimensions of the
  plot panels on the plot device will be proportional
  to the spatial dimensions of the
  corresponding components of \code{x}. This means that the
  objects will be plotted at \emph{approximately} equal physical scales.
  If these objects have very different spatial sizes,
  the plot command could fail (when it tries
  to plot the smaller objects at a tiny scale), with an error
  message that the figure margins are too large.

  The objects will be plotted at \emph{exactly} equal
  physical scales, and \emph{exactly} aligned on the device,
  under the following conditions:
  \itemize{
    \item
    every component of \code{x} is a spatial object
    whose position can be shifted by \code{\link{shift}};
    \item
    \code{panel.begin} and \code{panel.end} are either
    \code{NULL} or they are spatial objects 
    whose position can be shifted by \code{\link{shift}};
    \item
    \code{adorn.left}, 
    \code{adorn.right}, 
    \code{adorn.top} and 
    \code{adorn.bottom} are all \code{NULL}.
  }
}
\seealso{
  \code{\link{split.ppp}},
  \code{\link{plot.ppp}},
  \code{\link{ppp.object}},
  \code{\link{plot.listof}},
}
\section{Error messages}{
  If the error message \sQuote{Figure margins too large}
  occurs, ensure that \code{equal.scales=FALSE}
  and increase the values of \code{mar.panel}.
}
\examples{
# Multitype point pattern
 plot(split(amacrine))
 plot(split(amacrine), main="", 
     panel.begin=function(i, y, ...) { plot(density(y), ribbon=FALSE, ...) })
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{hplot}
