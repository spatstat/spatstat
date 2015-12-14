\name{rnoise}
\alias{rnoise}
\title{
  Random Pixel Noise
}
\description{
  Generate a pixel image whose pixel values are random numbers
  following a specified probability distribution.
}
\usage{
rnoise(rgen = runif, w = square(1), \dots)
}
\arguments{
  \item{rgen}{
    Random generator for the pixel values.
    A function in the \R language.
  }
  \item{w}{
    Window (region or pixel raster) in which to generate the image.
    Any data acceptable to \code{\link{as.mask}}.
  }
  \item{\dots}{
    Arguments, matched by name,
    to be passed to \code{rgen} to specify the parameters of the
    probability distribution, or passed to \code{\link{as.mask}}
    to control the pixel resolution.
  }
}
\details{
  The argument \code{w} could be a window (class \code{"owin"}),
  a pixel image (class \code{"im"}) or other data. It is 
  first converted to a binary mask by \code{\link{as.mask}}
  using any relevant arguments in \code{\dots}.
  
  Then each pixel inside the window (i.e. with logical value
  \code{TRUE} in the mask) is assigned a random numerical value
  by calling the function \code{rgen}.

  The function \code{rgen} would typically be one of the standard
  random variable generators like \code{\link{runif}} (uniformly
  distributed random values) or \code{\link{rnorm}} (Gaussian random
  values). Its first argument \code{n} is the number of values to be
  generated. Other arguments to \code{rgen} must be matched by name.
}
\value{
  A pixel image (object of class \code{"im"}).
}
\author{
Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}


Rolf Turner \email{r.turner@auckland.ac.nz}

and Ege Rubak \email{rubak@math.aau.dk}

}
\seealso{
  \code{\link{as.mask}},
  \code{\link{as.im}},
  \code{\link[stats]{Distributions}}.
}
\examples{
  plot(rnoise(), main="Uniform noise")
  plot(rnoise(rnorm, dimyx=32, mean=2, sd=1),
       main="White noise")
}
\keyword{spatial}
\keyword{datagen}

