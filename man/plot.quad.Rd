\name{plot.quad}
\alias{plot.quad}
\title{plot a Spatial Quadrature Scheme}
\description{
  Plot a two-dimensional spatial quadrature scheme.
}
\usage{
  \method{plot}{quad}(x, ..., main, dum=list())
}
\arguments{
  \item{x}{
    The spatial quadrature scheme to be plotted.
    An object of class \code{"quad"}.
  }
  \item{\dots}{
    extra arguments controlling the plotting
    of the data points of the quadrature scheme.
  }
  \item{dum}{
    list of extra arguments controlling the plotting of the dummy points of the
    quadrature scheme. See below.
  }
  \item{main}{
    text to be displayed as a title above the plot.
  }
}
\value{
  \code{NULL}.
}
\details{
  This is the \code{plot} method for 
  quadrature schemes (objects of class \code{"quad"},
  see \code{\link{quad.object}}).

  First the data points of the quadrature scheme
  are plotted (in their observation window) using 
  \code{\link{plot.ppp}} with any arguments specified in \code{...}

  Then the dummy points of the quadrature scheme are plotted
  using \code{\link{plot.ppp}} with any arguments specified in
  \code{dum}.

  By default the dummy points are superimposed onto the plot of data
  points. This can be overridden by including the argument
  \code{add=FALSE} in the list \code{dum} as shown in the examples.
  In this case the data and dummy point patterns are plotted separately.

  See \code{\link{par}} and \code{\link{plot.ppp}}
  for other possible arguments controlling the plots.
}
\seealso{
  \code{\link{quad.object}},
  \code{\link{plot.ppp}},
  \code{\link{par}}
}
\examples{
   data(nztrees)
   Q <- quadscheme(nztrees)

   plot(Q, main="NZ trees: quadrature scheme")

   oldpar <- par(mfrow=c(2,1))
   plot(Q, main="NZ trees", dum=list(add=FALSE))
   par(oldpar)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
 
 \keyword{hplot}
