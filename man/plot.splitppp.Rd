\name{plot.splitppp}
\alias{plot.splitppp}
\title{Plot a List of Point Patterns}
\description{
  Plots a list of point patterns.
}
\usage{
  \method{plot}{splitppp}(x, \dots, main, arrange=TRUE,
   nrows=NULL, ncols=NULL, main.panel=NULL, mar.panel=c(2,1,1,2),
   panel.begin=NULL, panel.end=NULL, panel.args=NULL, plotcommand="plot",
   adorn.left=NULL, adorn.right=NULL, adorn.top=NULL, adorn.bottom=NULL,
   adorn.size=0.2)
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
}
\seealso{
  \code{\link{split.ppp}},
  \code{\link{plot.ppp}},
  \code{\link{ppp.object}}
}
\examples{
# Multitype point pattern
 data(amacrine)
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
