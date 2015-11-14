\name{as.data.frame.im}
\alias{as.data.frame.im}
\title{Convert Pixel Image to Data Frame}
\description{
  Convert a pixel image to a data frame
}
\usage{
  \method{as.data.frame}{im}(x, ...)
}
\arguments{
  \item{x}{A pixel image (object of class \code{"im"}).}
  \item{\dots}{Further arguments passed to
    \code{\link[base:as.data.frame]{as.data.frame.default}}
    to determine the row names and other features.}
}
\details{
  This function takes the pixel image \code{x}
  and returns a data frame with three columns
  containing the pixel coordinates and the pixel values. 
}
\value{
  A data frame
}
\examples{
   # artificial image
   Z <- setcov(square(1))

   Y <- as.data.frame(Z)

   head(Y)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{methods}
