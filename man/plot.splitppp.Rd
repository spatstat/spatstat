\name{plot.splitppp}
\alias{plot.splitppp}
\title{Plot a List of Point Patterns}
\description{
  Plots a list of point patterns.
}
\usage{
  \method{plot}{splitppp}(x, \dots, main)
}
\arguments{
  \item{x}{
    A named list of point patterns,
    typically obtained from \code{\link{split.ppp}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.listof}} which control the
    layout of the plot panels, their appearance, and the
    plot behaviour in individual plot panels.
  }
  \item{main}{
    Optional main title for the plot.
  }
}
\value{
  Null.
}
\details{
  This is the \code{plot} method for the class \code{"splitppp"}.
  It is typically used to plot the result of the function
  \code{\link{split.ppp}}.

  The argument \code{x} should be a named list of point patterns
  (objects of class \code{"ppp"}, see \code{\link{ppp.object}}).
  Each of these point patterns will be plotted in turn
  using \code{\link{plot.ppp}}.

  Plotting is performed by \code{\link{plot.listof}}.
}
\seealso{
  \code{\link{plot.listof}} for arguments controlling the plot.
  
  \code{\link{split.ppp}},
  \code{\link{plot.ppp}},
  \code{\link{ppp.object}}.
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
