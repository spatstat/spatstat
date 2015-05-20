\name{density.splitppp}
\alias{density.splitppp}
\alias{density.ppplist}
\title{Kernel Smoothed Intensity of Split Point Pattern}
\description{
  Compute a kernel smoothed intensity function for each
  of the components of a split point pattern,
  or each of the point patterns in a list.
}
\usage{
  \method{density}{splitppp}(x, \dots, se=FALSE)

  \method{density}{ppplist}(x, \dots, se=FALSE)
}
\arguments{
  \item{x}{
    Split point pattern (object of class \code{"splitppp"}
    created by \code{\link{split.ppp}}) to be smoothed.
    Alternatively a list of point patterns,
    of class \code{"ppplist"}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{density.ppp}} to control
    the smoothing, pixel resolution, edge correction etc.
  }
  \item{se}{
    Logical value indicating whether to compute standard errors as well.
  }
}
\value{
  A list of pixel images (objects of class \code{"im"})
  which can be plotted or printed;
  or a list of numeric vectors giving the values at specified points.

  If \code{se=TRUE}, the result is a list with two elements named
  \code{estimate} and \code{SE}, each of the format described above.
}
\details{
  This is a method for the generic function \code{density}.

  The argument \code{x} should be a list of point patterns,
  and should belong to one of the classes 
  \code{"ppplist"} or \code{"splitppp"}.
  
  Typically \code{x} is obtained by applying
  the function \code{\link{split.ppp}} to a point pattern \code{y}
  by calling \code{split(y)}. This splits the points of \code{y} into several
  sub-patterns.
  
  A kernel estimate of the intensity function of each of the
  point patterns is computed using \code{\link{density.ppp}}.

  The return value is usually a list, each of whose entries is a
  pixel image (object of class \code{"im"}). The return value
  also belongs to the class \code{"solist"} and can be plotted
  or printed.

  If the argument \code{at="points"} is given, the result is a list
  of numeric vectors giving the intensity values at the data points.

  If \code{se=TRUE}, the result is a list with two elements named
  \code{estimate} and \code{SE}, each of the format described above.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{im.object}}
}
\examples{
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
