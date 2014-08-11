\name{pixellate}
\Rdversion{1.1}
\alias{pixellate}
\title{
  Convert Spatial Object to Pixel Image
}
\description{
  Convert a spatial object to 
  a pixel image by measuring the amount of stuff
  in each pixel.
}
\usage{
pixellate(x, ...)
}
\arguments{
  \item{x}{
    Spatial object to be converted.
    A point pattern (object of class \code{"ppp"}),
    a window (object of class \code{"owin"}),
    a line segment pattern (object of class \code{"psp"}),
    or some other suitable data.
}
  \item{\dots}{
    Arguments passed to methods.
}
}
\details{
  The function \code{pixellate} converts a geometrical object \code{x} into
  a pixel image, by measuring the \emph{amount} of \code{x} that is
  inside each pixel.
  
  If \code{x} is a point pattern, \code{pixellate(x)} counts the
  number of points of \code{x} falling in each pixel.
  If \code{x} is a window, \code{pixellate(x)} measures the area
  of intersection of each pixel with the window.

  The function \code{pixellate} is generic, with methods for
  point patterns (\code{\link{pixellate.ppp}}),
  windows (\code{\link{pixellate.owin}}),
  and 
  line segment patterns (\code{\link{pixellate.psp}}),
  See the separate documentation for these methods.

  The related function \code{\link{as.im}}
  also converts \code{x} into a pixel image,
  but typically measures only the presence or absence of
  \code{x} inside each pixel.
}
\value{
  A pixel image (object of class \code{"im"}).
}
\seealso{
  \code{\link{pixellate.ppp}},
  \code{\link{pixellate.owin}},
  \code{\link{pixellate.psp}},
  \code{\link{as.im}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
