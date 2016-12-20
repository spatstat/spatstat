\name{as.data.frame.owin}
\alias{as.data.frame.owin}
\title{Convert Window to Data Frame}
\description{
  Converts a window object to a data frame.
}
\usage{
\method{as.data.frame}{owin}(x, \dots, drop=TRUE)
}
\arguments{
  \item{x}{
    Window (object of class \code{"owin"}).
  }
  \item{\dots}{Further arguments passed to
    \code{\link[base:as.data.frame]{as.data.frame.default}}
    to determine the row names and other features.
  }
  \item{drop}{
    Logical value indicating whether to discard pixels that are
    outside the window, when \code{x} is a binary mask.
  }
}
\details{
  This function returns a data frame specifying the coordinates of the
  window.
  
  If \code{x} is a binary mask window,
  the result is a data frame with columns \code{x} and \code{y}
  containing the spatial coordinates of each \emph{pixel}.
  If \code{drop=TRUE} (the default), only pixels inside the window are retained.
  If \code{drop=FALSE}, all pixels are retained, and the data frame has
  an extra column \code{inside} containing the logical value of each pixel
  (\code{TRUE} for pixels inside the window, \code{FALSE} for outside).
  
  If \code{x} is a rectangle or a polygonal window,
  the result is a data frame with columns \code{x} and \code{y}
  containing the spatial coordinates of the \emph{vertices} of the
  window. If the boundary consists of several polygons, the data frame
  has additional columns \code{id}, identifying which polygon is being
  traced, and \code{sign}, indicating whether the polygon is an
  outer or inner boundary (\code{sign=1} and \code{sign=-1} respectively).
}
\value{
  A data frame with columns named \code{x} and \code{y},
  and possibly other columns.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{as.data.frame.im}}
}
\examples{
   as.data.frame(square(1))

   holey <- owin(poly=list(
                        list(x=c(0,10,0), y=c(0,0,10)),
                        list(x=c(2,2,4,4), y=c(2,4,4,2))))
   as.data.frame(holey)
}
\keyword{spatial}
\keyword{methods}

