\name{plot.lppm}
\alias{plot.lppm}
\title{
  Plot a Fitted Point Process Model on a Linear Network
}
\description{
  Plots the fitted intensity of a point process model
  on a linear network. 
}
\usage{
  \method{plot}{lppm}(x, ..., type="trend")
}
\arguments{
  \item{x}{
    An object of class \code{"lppm"} representing a fitted point process
    model on a linear network.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.linim}} to control the
    plot.
  }
  \item{type}{
    Character string (either \code{"trend"} or \code{"cif"})
    determining whether to plot the fitted first order trend
    or the conditional intensity.
  }
}
\details{
  This function is the plot method for the class 
  \code{"lppm"}. It computes the fitted intensity of the
  point process model, and displays it using \code{\link{plot.linim}}.

  The default is to display intensity values as colours. Alternatively
  if the argument \code{style="width"} is given, intensity values are
  displayed as the widths of thick lines drawn over the network.
}
\value{
  Null.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
}
\seealso{
  \code{\link{lppm}},
  \code{\link{plot.linim}},
  \code{\link{methods.lppm}},
  \code{\link{predict.lppm}}.
}
\examples{
  X <- runiflpp(10, simplenet)  
  fit <- lppm(X ~x)
  plot(fit)
  plot(fit, style="width")
}
\keyword{spatial}
\keyword{models}
