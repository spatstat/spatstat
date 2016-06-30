\name{is.im}
\alias{is.im}
\title{Test Whether An Object Is A Pixel Image}
\description{
  Tests whether its argument is a pixel image
  (object of class \code{"im"}).
}
\usage{
is.im(x)
}
\arguments{
  \item{x}{Any object.}
}
\details{
  This function tests whether the argument \code{x} is a
  pixel image object of class \code{"im"}. For details of this
  class, see \code{\link{im.object}}.
  
  The object is determined to be an image if it inherits from
  class \code{"im"}.
}
\value{
  \code{TRUE} if \code{x} is a pixel image, otherwise \code{FALSE}.
}

\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{manip}
