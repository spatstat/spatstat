\name{anyNA.im}
\alias{anyNA.im}
\title{
  Check Whether Image Contains NA Values
}
\description{
  Checks whether any pixel values in a pixel image are \code{NA}
  (meaning that the pixel lies outside the domain of definition of the
  image). 
}
\usage{
  \method{anyNA}{im}(x, recursive = FALSE)
}
\arguments{
  \item{x}{
    A pixel image (object of class \code{"im"}).
  }
  \item{recursive}{
    Ignored.
  }
}
\details{
  The function \code{\link{anyNA}} is generic: \code{anyNA(x)} is a
  faster alternative to \code{any(is.na(x))}.

  This function \code{anyNA.im} is a method for the generic \code{anyNA}
  defined for pixel images. It returns the value \code{TRUE} if any of the pixel
  values in \code{x} are \code{NA}, and
  and otherwise returns \code{FALSE}.
}
\value{
  A single logical value.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{im.object}}
}
\examples{
  anyNA(as.im(letterR))
}
\keyword{spatial}
\keyword{methods}


