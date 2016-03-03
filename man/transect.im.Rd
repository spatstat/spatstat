\name{transect.im}
\alias{transect.im}
\title{
  Pixel Values Along a Transect
}
\description{
  Extract the pixel values of a pixel image at each point along
  a linear transect.
}
\usage{
transect.im(X, ..., from="bottomleft", to="topright",
             click=FALSE, add=FALSE)
}
\arguments{
  \item{X}{
    A pixel image (object of class \code{"im"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{from,to}{
    Optional. Start point and end point of the transect.
    Pairs of \eqn{(x,y)} coordinates
    in a format acceptable to \code{\link{xy.coords}},
    or keywords \code{"bottom"}, \code{"left"}, \code{"top"},
    \code{"right"}, \code{"bottomleft"} etc.
  }
  \item{click}{
    Optional.
    Logical value.
    If \code{TRUE}, the linear transect is determined interactively
    by the user, who clicks two points on the current plot.
  }
  \item{add}{
    Logical. If \code{click=TRUE}, this argument determines
    whether to perform interactive tasks on the current plot (\code{add=TRUE})
    or to start by plotting \code{X} (\code{add=FALSE}).
  }
}
\details{
  The pixel values of the image \code{X} along a line segment 
  will be extracted. The result is a function table (\code{"fv"} object)
  which can be plotted directly.

  If \code{click=TRUE}, then the user is prompted to click two points on
  the plot of \code{X}. These endpoints define the transect.

  Otherwise, the transect is defined by the endpoints
  \code{from} and \code{to}. The default is a diagonal transect from
  bottom left to top right of the frame.
}
\value{
  An object of class \code{"fv"} which can be plotted.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\seealso{
  \code{\link{im}}
}
\examples{
  Z <- density(redwood)
  plot(transect.im(Z))
  \dontrun{
   if(FALSE) {
    plot(transect.im(Z, click=TRUE))
   }
  }
}
\keyword{spatial}
\keyword{manip}
\keyword{iplot}
