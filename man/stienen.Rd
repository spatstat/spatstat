\name{stienen}
\alias{stienen}
\alias{stienenSet}
\title{
  Stienen Diagram
}
\description{
  Draw the Stienen diagram of a point pattern,
  or compute the region covered by the Stienen diagram.
}
\usage{
stienen(X, \dots, bg = "grey", border = list(bg = NULL))
stienenSet(X, edge=TRUE)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.ppp}}
    to control the plot.
  }
  \item{bg}{
    Fill colour for circles.
  }
  \item{border}{
    Either a list of arguments passed to \code{\link{plot.ppp}}
    to control the display of circles at the border of the diagram,
    or the value \code{FALSE}
    indicating that the border circles should not be plotted.
  }
  \item{edge}{
    Logical value indicating whether to include the circles at the
    border of the diagram.
  }
}
\details{
  The Stienen diagram of a point pattern (Stienen, 1982)
  is formed by drawing a circle around each point of the pattern,
  with diameter equal to the nearest-neighbour distance for that point.
  These circles do not overlap. If two points are nearest neighbours
  of each other, then the corresponding circles touch.

  \code{stienenSet(X)} computes the union of these circles and
  returns it as a window (object of class \code{"owin"}).

  \code{stienen(X)} generates a plot of the Stienen diagram of
  the point pattern \code{X}. By default, circles are shaded in grey
  if they lie inside the window of \code{X}, and are not shaded
  otherwise.
}
\value{
  The plotting function \code{stienen} returns \code{NULL}.

  The return value of \code{stienenSet} is a window (object of class
  \code{"owin"}).
}
\references{
  Stienen, H. (1982)
  \emph{Die Vergroeberung von Karbiden in reinen Eisen-Kohlenstoff
    Staehlen}.
  Dissertation, RWTH Aachen.
}
\seealso{
  \code{\link{nndist}},
  \code{\link{plot.ppp}}
}
\examples{
  Y <- stienenSet(cells)
  stienen(redwood)
  stienen(redwood, border=list(bg=NULL, lwd=2, cols="red"))
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\keyword{spatial}
\keyword{math}
\keyword{manip}
