\name{is.linim}
\alias{is.linim}
\title{Test Whether an Object is a Pixel Image on a Linear Network}
\description{
  Tests whether its argument is a pixel image on a linear network
  (object of class \code{"linim"}).
}
\usage{
is.linim(x)
}
\arguments{
  \item{x}{Any object.}
}
\details{
  This function tests whether the argument \code{x} is a
  pixel image on a linear network
  (object of class \code{"linim"}). 
  
  The object is determined to be an image if it inherits from
  class \code{"linim"}.
}
\value{
  \code{TRUE} if \code{x} is a pixel image on a linear network,
  otherwise \code{FALSE}.
}

\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{manip}
