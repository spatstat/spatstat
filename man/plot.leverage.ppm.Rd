\name{plot.leverage.ppm}
\alias{plot.leverage.ppm}
\title{
  Plot Leverage Function
}
\description{
  Plots a leverage function that has been computed by \code{\link{leverage.ppm}}.
}
\usage{
 \method{plot}{leverage.ppm}(x, ..., showcut=TRUE, col.cut=par("fg"))
}
\arguments{
  \item{x}{
    Leverage measure (object of class \code{"leverage.ppm"}) computed by
    \code{\link{leverage.ppm}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.im}}
    or \code{\link{contour.im}} controlling the plot.
  }
  \item{showcut}{
    Logical. If \code{TRUE}, a contour line is plotted at the
    level equal to the theoretical mean of the leverage.
  }
  \item{col.cut}{
    Optional colour for the contour line.
  }
}
\details{
  This is the plot method for objects of class \code{"leverage.ppm"}.
  These objects are computed by the command \code{\link{leverage.ppm}}.

  Use the argument \code{clipwin} to restrict the plot to a subset
  of the full data.
}
\value{
  None.
}
\references{
  Baddeley, A., Chang, Y.M. and Song, Y. (2013)
  Leverage and influence diagnostics for spatial point process models.
  \emph{Scandinavian Journal of Statistics} \bold{40}, 86--104.
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

