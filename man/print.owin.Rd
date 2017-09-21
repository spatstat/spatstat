\name{print.owin}
\alias{print.owin}
\title{Print Brief Details of a Spatial Window}
\description{
  Prints a very brief description of a window object.
}
\usage{
  \method{print}{owin}(x, \dots, prefix="window: ")
}
\arguments{
  \item{x}{Window (object of class \code{"owin"}).}
  \item{\dots}{Ignored.}
  \item{prefix}{Character string to be printed at the start of the output.}
}
\details{
  A very brief description of the window \code{x} is printed.

  This is a method for the generic function \code{\link{print}}.
}
\seealso{
  \code{\link{print}},
  \code{\link{print.ppp}},
  \code{\link{summary.owin}}
}
\examples{
  owin()  # the unit square

  data(demopat)
  W <- Window(demopat)
  W                    # just says it is polygonal
  as.mask(W)           # just says it is a binary image 

}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{print}

