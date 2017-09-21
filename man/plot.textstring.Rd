\name{plot.textstring}
\alias{plot.textstring}
\title{Plot a Text String}
\description{Plots an object of class \code{"textstring"}.}
\usage{
\method{plot}{textstring}(x, \dots, do.plot = TRUE)
}
\arguments{
  \item{x}{
    Object of class \code{"textstring"} to be plotted.
    This object is created by the command \code{\link{textstring}}.
  }
  \item{\dots}{
    Additional graphics arguments passed to
    \code{\link[graphics]{text}} to control the plotting of text.
  }
  \item{do.plot}{
    Logical value indicating whether to actually plot the text.
  }
}
\details{
  The argument \code{x} should be an object of class \code{"textstring"} 
  created by the command \code{\link{textstring}}. 

  This function displays the text using 
  \code{\link[graphics]{text}}.
}
\value{
  A window (class \code{"owin"}) enclosing the plotted graphics.
}
\examples{
  W <- Window(humberside)
  te <- textstring(centroid.owin(W), txt="Humberside", cex=2.5)
  plot(layered(W, te), main="")
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\seealso{
  \code{\link{onearrow}},
  \code{\link{yardstick}}
}
\keyword{spatial}
\keyword{hplot}
