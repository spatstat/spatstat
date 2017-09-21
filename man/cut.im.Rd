\name{cut.im}
\alias{cut.im}
\title{Convert Pixel Image from Numeric to Factor}
\description{
  Transform the values of a pixel image
  from numeric values into a factor.
}
\usage{
  \method{cut}{im}(x, \dots)
}
\arguments{
  \item{x}{
    A pixel image.
    An object of class \code{"im"}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{cut.default}}.
    They determine the breakpoints for the mapping from numerical values to
    factor values. See \code{\link{cut.default}}.
  }
} 
\value{
  A pixel image (object of class \code{"im"}) with 
  pixel values that are a factor.
  See \code{\link{im.object}}.
}
\details{
  This simple function applies the generic \code{\link{cut}} operation
  to the pixel values of the image \code{x}. The range of pixel values
  is divided into several intervals, and each
  interval is associated with a level of a factor. 
  The result is another pixel image,
  with the same window and pixel grid as
  \code{x}, but with the numeric value of each pixel discretised
  by replacing it by the factor level. 

  This function is a convenient
  way to inspect an image and to obtain summary statistics.
  See the examples.

  To select a subset of an image, use the subset operator
  \code{\link{[.im}} instead.
}
\seealso{
  \code{\link{cut}},
  \code{\link{im.object}}
}
\examples{
  # artificial image data
  Z <- setcov(square(1))

  Y <- cut(Z, 3)
  Y <- cut(Z, breaks=seq(0,1,length=5))

  # cut at the quartiles
  # (divides the image into 4 equal areas)
  Y <- cut(Z, quantile(Z))
}

\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{methods}
