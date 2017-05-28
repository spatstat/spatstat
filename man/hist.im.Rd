\name{hist.im}
\alias{hist.im}
\title{Histogram of Pixel Values in an Image}
\description{
  Computes and displays a histogram of the pixel values in a pixel image.
  The \code{hist} method for class \code{"im"}.
}
\usage{
  \method{hist}{im}(x, \dots, probability=FALSE, xname)
}
\arguments{
  \item{x}{A pixel image (object of class \code{"im"}).}
  \item{\dots}{Arguments passed to \code{\link{hist.default}}
    or \code{\link{barplot}}.}
  \item{probability}{Logical. If \code{TRUE}, the histogram will be
    normalised to give probabilities or probability densities.
  }
  \item{xname}{Optional. Character string to be used as the
  name of the dataset \code{x}.
  }
}
\details{
  This function computes and (by default) displays a histogram
  of the pixel values in the image \code{x}.

  An object of class \code{"im"}
  describes a pixel image. See \code{\link{im.object}})
  for details of this class.

  The function \code{hist.im} is a method for the generic
  function \code{\link{hist}} for the class \code{"im"}. 
  
  Any arguments in \code{...} are passed to \code{\link{hist.default}}
  (for numeric valued images) or \code{\link{barplot}} (for factor or
  logical images).
  For example, such arguments control the axes, and may be used to
  suppress the plotting.
}
\value{
  For numeric-valued images, an object of class \code{"histogram"} as returned
  by \code{\link[graphics:hist]{hist.default}}. This object can be
  plotted.
  
  For factor-valued or logical images, an object of class
  \code{"barplotdata"}, which can be plotted.
  This is a list with components
  called \code{counts} (contingency table of counts of the numbers of
  pixels taking each possible value), \code{probs} (corresponding relative
  frequencies) and \code{mids} (graphical \eqn{x}-coordinates of the
  midpoints of the bars in the barplot). 
}
\seealso{
  \code{\link{spatialcdf}} for the cumulative distribution function
  of an image.

  \code{\link{hist}},
  \code{\link{hist.default}},
  \code{\link{barplot}}.

  For other statistical graphics such as Q-Q plots,
  use \code{X[]} to extract the pixel values of image \code{X},
  and apply the usual statistical graphics commands.
  
  For information about pixel images see
  \code{\link{im.object}},
  \code{\link{summary.im}}.
}
\examples{
  X <- as.im(function(x,y) {x^2}, unit.square())
  hist(X)
  hist(cut(X,3))
}
\author{
  \adrian
  and \rolf
}
\keyword{spatial}
\keyword{methods}
