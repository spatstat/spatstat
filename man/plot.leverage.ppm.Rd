\name{plot.leverage.ppm}
\alias{plot.leverage.ppm}
\title{
  Plot Leverage Function
}
\description{
  Plots a leverage function that has been computed by \code{\link{leverage.ppm}}.
}
\usage{
 \method{plot}{leverage.ppm}(x, ..., showcut=TRUE)
}
\arguments{
  \item{x}{
    Leverage measure (object of class \code{"leverage.ppm"}) computed by
    \code{\link{leverage.ppm}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat]{plot.im}}
    controlling the image plot.
  }
  \item{showcut}{
    Logical. If \code{TRUE}, a contour line is plotted at the
    level equal to the theoretical mean of the leverage.
  }
}
\details{
  This is the plot method for objects of class \code{"leverage.ppm"}.
  These objects are computed by the command \code{\link{leverage.ppm}}.
}
\value{
  None.
}
\references{
  Baddeley, A. and Chang, Y.M. and Song, Y. (2011)
  Leverage and influence diagnostics for spatial point process models.
  \emph{Scandinavian Journal of Statistics}, in press.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{leverage.ppm}}
}
\examples{
   X <- rpoispp(function(x,y) { exp(3+3*x) })
   fit <- ppm(X, ~x+y)
   plot(leverage(fit))
}
\keyword{spatial}
\keyword{models}

