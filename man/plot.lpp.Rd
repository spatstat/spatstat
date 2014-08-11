\name{plot.lpp}
\alias{plot.lpp}
\title{
  Plot Point Pattern on Linear Network
}
\description{
  Plots a point pattern on a linear network.
  Plot method for the class \code{"lpp"} of point patterns on a linear network.
}
\usage{
\method{plot}{lpp}(x, \dots, main, add = FALSE,
                   use.marks=TRUE, which.marks=NULL,
                   show.all = !add, do.plot = TRUE, multiplot=TRUE) 
}
\arguments{
  \item{x}{
    Point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{plot.linnet}}
    or \code{\link{plot.ppp}}.
  }
  \item{main}{
    Main title for plot.
  }
  \item{add}{
    Logical value indicating whether the plot is to be added to the
    existing plot (\code{add=TRUE}) or whether a new plot should be
    initialised (\code{add=FALSE}, the default).
  }
  \item{use.marks}{
    logical flag; if \code{TRUE}, plot points using a different
    plotting symbol for each mark;
    if \code{FALSE}, only the locations of the points will be plotted,
    using \code{\link{points}()}.
  }
  \item{which.marks}{
    Index determining which column of marks to use,
    if the marks of \code{x} are a data frame.
    A character or integer vector identifying one or more
    columns of marks. 
    If \code{add=FALSE} then
    the default is to plot all columns of marks, in a series of separate
    plots.
    If \code{add=TRUE} then only one column of marks can be plotted,
    and the default is \code{which.marks=1}
    indicating the first column of marks.
  }
  \item{show.all}{
    Logical value indicating whether to plot everything
    including the main title and the observation window of \code{x}.
  }
  \item{do.plot}{
    Logical value determining whether to actually perform the plotting.
  }
  \item{multiplot}{
    Logical value giving permission to display multiple plots.
  }
}
\details{
  The linear network is plotted by \code{\link{plot.linnet}},
  then the points are plotted by \code{\link{plot.ppp}}.

  Commonly-used arguments include:
  \itemize{
    \item \code{col} and \code{lwd} for the colour and width of lines
    in the linear network
    \item \code{cols} for the colour or colours of the points
    \item \code{chars} for the plot characters representing different
    types of points
    \item \code{legend} and \code{leg.side} to control the graphics
    legend
    \item \code{show.window=FALSE} to suppress the box drawn around
    the linear network.
  }
}
\value{
  (Invisible) object of class \code{"symbolmap"}
  giving the correspondence between 
  mark values and plotting characters.
}
\seealso{
  \code{\link{lpp}}.
  
  See \code{\link{plot.ppp}} for options for representing the points.
}
\examples{
  plot(chicago, cols=1:6)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{hplot}
