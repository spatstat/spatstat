\name{hist.funxy}
\alias{hist.funxy}
\title{Histogram of Values of a Spatial Function}
\description{
  Computes and displays a histogram of the values of a spatial function
  of class \code{"funxy"}.
}
\usage{
  \method{hist}{funxy}(x, \dots, xname)
}
\arguments{
  \item{x}{A pixel image (object of class \code{"funxy"}).}
  \item{\dots}{
    Arguments passed to \code{\link{as.im}}
    or \code{\link{hist.im}}.
  }
  \item{xname}{
    Optional. Character string to be used as the
    name of the dataset \code{x}.
  }
}
\details{
  This function computes and (by default) displays a histogram
  of the values of the function \code{x}.

  An object of class \code{"funxy"} 
  describes a function of spatial location. It is a \code{function(x,y,..)}
  in the \R language, with additional attributes.

  The function \code{hist.funxy} is a method for the generic
  function \code{\link{hist}} for the class \code{"funxy"}.

  The function is first converted to a pixel image using \code{\link{as.im}},
  then \code{\link{hist.im}} is called to produce the histogram.

  Any arguments in \code{...} are passed to \code{\link{as.im}}
  to determine the pixel resolution, 
  or to \code{\link{hist.im}} to determine the histogram breaks
  and to control or suppress plotting.
  Useful arguments include \code{W} for the spatial domain,
  \code{eps,dimyx} for pixel resolution, \code{main} for the main title.
}
\value{
  An object of class \code{"histogram"} as returned
  by \code{\link[graphics:hist]{hist.default}}. This object can be
  plotted.
}
\seealso{
  \code{\link{spatialcdf}} for the cumulative distribution function
  of an image or function.

  \code{\link{hist}},
  \code{\link{hist.default}}.

  For other statistical graphics such as Q-Q plots,
  use \code{as.im(X)[]} to extract the pixel values of image \code{X},
  and apply the usual statistical graphics commands.
}
\examples{
  f <- funxy(function(x,y) {x^2}, unit.square())
  hist(f)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
