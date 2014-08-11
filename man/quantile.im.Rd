\name{quantile.im}
\alias{quantile.im}
\title{Sample Quantiles of Pixel Image}
\description{
Compute the sample quantiles of the pixel values of a
given pixel image.
}
\usage{
\method{quantile}{im}(x, \dots)
}
\arguments{
\item{x}{
A pixel image.
An object of class \code{"im"}.
}
\item{\dots}{
Optional arguments passed to \code{\link{quantile.default}}.
They determine the probabilities for which quantiles should be
computed. See \code{\link{quantile.default}}.
}
} 
\value{
A vector of quantiles. 
}
\details{
This simple function applies the generic \code{\link{quantile}} operation
to the pixel values of the image \code{x}. 

This function is a convenient
way to inspect an image and to obtain summary statistics.
See the examples.
}
\seealso{
\code{\link{quantile}},
\code{\link{cut.im}},
\code{\link{im.object}}
}
\examples{
# artificial image data
Z <- setcov(square(1))

# find the quartiles
quantile(Z)

# find the deciles
quantile(Z, probs=(0:10)/10)
}

\author{Adrian Baddeley
\email{Adrian.Baddeley@csiro.au}
\url{http://www.maths.uwa.edu.au/~adrian/}
and Rolf Turner
\email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{univar}
