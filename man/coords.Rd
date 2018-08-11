\name{coords}
\Rdversion{1.1}
\alias{coords}
\alias{coords.ppp}
\alias{coords.ppx}
\alias{coords<-}
\alias{coords<-.ppp}
\alias{coords<-.ppx}
\title{
  Extract or Change Coordinates of a Spatial or Spatiotemporal Point Pattern
}
\description{
  Given any kind of spatial or space-time point pattern,
  this function extracts the
  (space and/or time and/or local) coordinates of the points
  and returns them as a data frame.
}
\usage{
  coords(x, ...)
  \method{coords}{ppp}(x, ...)
  \method{coords}{ppx}(x, ..., spatial = TRUE, temporal = TRUE, local=TRUE)
  coords(x, ...) <- value
  \method{coords}{ppp}(x, ...) <- value
  \method{coords}{ppx}(x, ..., spatial = TRUE, temporal = TRUE, local=TRUE) <- value
}
\arguments{
  \item{x}{
    A point pattern: either a two-dimensional point pattern
    (object of class \code{"ppp"}), a three-dimensional point pattern
    (object of class \code{"pp3"}), or a 
    general multidimensional space-time point pattern 
    (object of class \code{"ppx"}).
  }
  \item{\dots}{
    Further arguments passed to methods.
  }
  \item{spatial,temporal,local}{
    Logical values indicating whether to extract spatial, temporal
    and local coordinates, respectively. The default is to return all
    such coordinates. (Only relevant to \code{ppx} objects).
  }
  \item{value}{
    New values of the coordinates. A numeric vector with one entry for each
    point in \code{x}, or a numeric matrix or data frame with one row
    for each point in \code{x}.
  }
}
\details{
  The function \code{coords} extracts the coordinates from
  a point pattern. The function \code{coords<-} replaces the coordinates
  of the point pattern with new values.
  
  Both functions \code{coords} and \code{coords<-} are generic, with methods for
  the classes \code{"ppp"}) and \code{"ppx"}.
  An object of class \code{"pp3"} also inherits from \code{"ppx"} and
  is handled by the method for \code{"ppx"}.
}
\value{
  \code{coords} returns a \code{data.frame} with one row for each point,
  containing the coordinates.
  \code{coords<-} returns the altered point pattern.
}
\author{\adrian
  
  
  and \rolf
  
}
\seealso{
  \code{\link{ppx}},
  \code{\link{pp3}},
  \code{\link{ppp}},
  \code{as.hyperframe.ppx},
  \code{as.data.frame.ppx}.
}
\examples{
   df <- data.frame(x=runif(4),y=runif(4),t=runif(4))
   X <- ppx(data=df, coord.type=c("s","s","t"))
   coords(X)
   coords(X, temporal=FALSE)
   coords(X) <- matrix(runif(12), ncol=3)
}
\keyword{spatial}
\keyword{manip}
