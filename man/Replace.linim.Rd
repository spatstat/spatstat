\name{Replace.linim}
\alias{[<-.linim}
\title{Reset Values in Subset of Image on Linear Network}
\description{
  Reset the values in a subset of a pixel image on a linear network.
}
\usage{
  \method{[}{linim}(x, i, j) <- value
}
\arguments{
  \item{x}{
    A pixel image on a linear network.
    An object of class \code{"linim"}.
  }
  \item{i}{
    Object defining the subregion or subset to be replaced.
    Either a spatial window (an object of class \code{"owin"}), or a
    pixel image with logical values, or a point pattern (an object
    of class \code{"ppp"}), or any type of index that applies to a
    matrix, or something that can be converted to a point pattern
    by \code{\link{as.ppp}} (using the window of \code{x}).
  }
  \item{j}{
    An integer or logical vector serving as the column index
    if matrix indexing is being used.  Ignored if \code{i} is
    appropriate to some sort of replacement \emph{other than}
    matrix indexing.
  }
  \item{value}{
     Vector, matrix, factor or pixel image
     containing the replacement values.
     Short vectors will be recycled.
  }
} 
\value{
  The image \code{x} with the values replaced.
}
\details{
  This function changes some of the pixel values in a
  pixel image. The image \code{x} must be an object of class
  \code{"linim"} representing a pixel image on a linear network.

  The pixel values are replaced according to the rules
  described in the help for \code{\link{[<-.im}}.
  Then the auxiliary data are updated.
}
\seealso{
  \code{\link{[<-.im}}.
}
\examples{
  # make a function
  Y <- as.linim(distfun(runiflpp(5, simplenet)))
  # replace some values
  B <- square(c(0.25, 0.55))
  Y[B] <- 2
  plot(Y, main="")
  plot(B, add=TRUE, lty=3)
  X <- runiflpp(4, simplenet)
  Y[X] <- 5
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
