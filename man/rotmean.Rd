\name{rotmean}
\alias{rotmean}
\title{
  Rotational Average of a Pixel Image
}
\description{
  Compute the average pixel value over all rotations of the image
  about the origin, as a function of distance from the origin.
}
\usage{
rotmean(X, ..., origin, Xname)
}
\arguments{
  \item{X}{
    A pixel image.
  }
  \item{\dots}{
    Ignored.
  }
  \item{origin}{
    Optional. Origin about which the rotations should be performed.
    Either a numeric vector or a character string as described
    in the help for \code{\link{shift.owin}}.
  }
  \item{Xname}{
    Optional name for \code{X} to be used in the
    function labels.
  }
}
\details{
  This command computes, for each possible distance \eqn{r},
  the average pixel value of the pixels lying at
  distance \eqn{r} from the origin. Kernel smoothing is used
  to obtain a smooth function of \eqn{r}.
}
\value{
  An object of class \code{"fv"} giving the mean pixel value of \code{X}
  as a function of distance from the origin.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\examples{
  if(interactive()) {
    Z <- setcov(square(1))
    plot(rotmean(Z))
  } else {
    Z <- setcov(square(1), dimyx=32)
    f <- rotmean(Z)
  }
}
\keyword{spatial}
\keyword{math}
