\name{plot.linnet}
\alias{plot.linnet}
\title{
  Plot a linear network
}
\description{
  Plots a linear network
}
\usage{
 \method{plot}{linnet}(x, ..., main=NULL, add=FALSE, 
                               vertices=FALSE, window=FALSE)
}
\arguments{
  \item{x}{
    Linear network (object of class \code{"linnet"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.psp}}
    controlling the plot.
  }
  \item{main}{
    Main title for plot. Use \code{main=""} to suppress it.
  }
  \item{add}{
    Logical. If code{TRUE}, superimpose the graphics
    over the current plot. If \code{FALSE}, generate a new plot.
  }
  \item{vertices}{
    Logical. Whether to plot the vertices as well.
  }
  \item{window}{
    Logical. Whether to plot the window containing the linear network.
  }
}
\details{
  This is the plot method for class \code{"linnet"}.
}
\value{
  Null.
}
\author{
  Ang Qi Wei \email{aqw07398@hotmail.com} and
  Adrian Baddeley \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
   \code{\link{linnet}}
}
\examples{
   example(linnet)
   plot(letterA)
}
\keyword{spatial}
