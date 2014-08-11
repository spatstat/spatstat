\name{density.splitppp}
\alias{density.splitppp}
\title{Kernel Smoothed Intensity of Split Point Pattern}
\description{
  Compute a kernel smoothed intensity function for each
  of the components of a split point pattern.
}
\usage{
  \method{density}{splitppp}(x, \dots)
}
\arguments{
  \item{x}{
    Split point pattern (object of class \code{"splitppp"}
    created by \code{\link{split.ppp}}) to be smoothed.
  }
  \item{\dots}{
    Arguments passed to \code{\link{density.ppp}} to control
    the smoothing, pixel resolution, edge correction etc.
  }
}
\value{
  A list of pixel images (objects of class \code{"im"}).
  Can be plotted or printed.
}
\details{
  This is a method for the generic function \code{density}.

  The argument \code{x} should be an object of class
  \code{"splitppp"}, effectively a list of point patterns.
  
  Typically \code{x} is obtained by applying
  the function \code{\link{split.ppp}} to a point pattern \code{y}
  by calling \code{split(y)}. This splits the points of \code{y} into several
  sub-patterns.
  
  A kernel estimate of the intensity function of each of the
  point patterns is computed using \code{\link{density.ppp}}.

  The return value is a list, each of whose entries is a
  pixel image (object of class \code{"im"}). The return value
  also belongs to the class \code{"listof"} and can be plotted
  or printed.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{im.object}}
}
\examples{
  data(amacrine)
  Z <- density(split(amacrine), 0.05)
  plot(Z)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
