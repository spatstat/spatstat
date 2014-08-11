\name{plot.listof}
\alias{plot.listof}
\title{Plot a List of Things}
\description{
  Plots a list of things
}
\usage{
  \method{plot}{listof}(x, \dots, main, arrange=TRUE,
   nrows=NULL, ncols=NULL, main.panel=NULL, mar.panel=c(2,1,1,2),
   panel.begin=NULL, panel.end=NULL, panel.args=NULL, plotcommand="plot",
   adorn.left=NULL, adorn.right=NULL, adorn.top=NULL, adorn.bottom=NULL,
   adorn.size=0.2)
}
\arguments{
  \item{x}{
    An object of the class \code{"listof"}.
    Essentially a list of objects.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot}} when generating each
    plot panel.
  }
  \item{main}{
    Overall heading for the plot.
  }
  \item{arrange}{
    Logical flag indicating whether to plot the objects
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
    giving the headings for each of the objects.
  }
  \item{mar.panel}{
    Value of the graphics parameter \code{mar} controlling the size
    of the margins outside each plot panel.
    See the help file for \code{\link{par}}.
  }
  \item{panel.begin,panel.end}{
    Optional. Functions that will be executed before and after each panel is
    plotted. See Details.
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
}
\value{
  Null.
}
\details{
  This is the \code{plot} method for the class \code{"listof"}.

  An object of class \code{"listof"} (defined in the base R package) represents
  a list of objects, all belonging to a common class.
  The base R package defines a method for printing these objects,
  \code{\link{print.listof}}, but does not define a method for \code{plot}.
  So here we have provided a method for \code{plot}.

  In the \pkg{spatstat} package, the function
  \code{\link{density.splitppp}} produces an object of
  class \code{"listof"}, essentially a list of pixel images.
  These images can be plotted in a nice arrangement
  using \code{plot.listof}. See the Example.

  The arguments \code{panel.begin} and \code{panel.end}
  may be functions that will be executed before and after each panel is plotted.
  They will be called as \code{panel.begin(i, y, main=main.panel[i])}
  and \code{panel.end(i, y, add=TRUE)}.

  Alternatively, \code{panel.begin} and \code{panel.end} may be objects
  of some class that can be plotted
  with the generic \code{plot} command. They will be plotted before and
  after each panel is plotted.

  If all entries of \code{x} are pixel images,
  the function \code{\link{image.listof}} is called to control
  the plotting. The arguments \code{equal.ribbon} and \code{col}
  can be used to determine the colour map or maps applied.
}
\seealso{
  \code{\link{print.listof}},
  \code{\link{contour.listof}},
  \code{\link{image.listof}},
  \code{\link{density.splitppp}}
}
\examples{
# Multitype point pattern
 data(amacrine)
 plot(D <- density(split(amacrine)))
 plot(D, main="", equal.ribbon=TRUE,
      panel.end=function(i,y,...){contour(y, ...)})
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{hplot}
