\name{harmonise.owin}
\alias{harmonise.owin}
\alias{harmonize.owin}
\title{Make Windows Compatible}
\description{
  Convert several windows to a common pixel raster.
}
\usage{
\method{harmonise}{owin}(\dots)

\method{harmonize}{owin}(\dots)
}
\arguments{
  \item{\dots}{
    Any number of windows (objects of class \code{"owin"})
    or data which can be converted to windows by \code{\link{as.owin}}.
  } 
}
\details{
  This function makes any number of windows compatible,
  by converting them all to a common pixel grid.

  This only has an effect if one of the windows is a binary mask.
  If all the windows are rectangular or polygonal, they are
  returned unchanged.

  The command \code{\link{harmonise}} is generic. This is the
  method for objects of class \code{"owin"}.
  
  Each argument must be a window (object of class \code{"owin"}),
  or data that can be converted to a window by \code{\link{as.owin}}.
  
  The common pixel grid is determined by inspecting all the windows
  in the argument list, computing the bounding box of all the
  windows, then finding the binary mask with the finest spatial resolution, 
  and extending its pixel grid to cover the bounding box. 

  The return value is a list with entries corresponding to the input
  arguments.
  If the arguments were named (\code{name=value}) then the return value
  also carries these names.

  If you just want to determine the appropriate pixel resolution,
  without converting the windows, use \code{\link{commonGrid}}.
}
\value{
  A list of windows, of length equal to the number of arguments
  \code{\dots}. The list belongs to the class \code{"solist"}.
}
\author{\adrian
  
  
  and \rolf
  
}
\examples{
   harmonise(X=letterR,
             Y=grow.rectangle(Frame(letterR), 0.2),
             Z=as.mask(letterR, eps=0.1),
             V=as.mask(letterR, eps=0.07))
}
\seealso{
  \code{\link{commonGrid}},
  \code{\link{harmonise.im}},
  \code{\link{as.owin}}
}
\keyword{spatial}
\keyword{manip}
