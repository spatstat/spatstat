\name{plot.rppm}
\alias{plot.rppm}
\title{
  Plot a Recursively Partitioned Point Process Model
}
\description{
  Given a model which has been fitted to point pattern data
  by recursive partitioning, plot the partition tree
  or the fitted intensity.
}
\usage{
\method{plot}{rppm}(x, \dots, what = c("tree", "spatial"))
}
\arguments{
  \item{x}{
    Fitted point process model of class \code{"rppm"}
    produced by the function \code{\link{rppm}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link[rpart]{plot.rpart}}
    or \code{\link{plot.im}}
    controlling the appearance of the plot.
  }
  \item{what}{
    Character string (partially matched) specifying whether to plot the
    partition tree or the fitted intensity.
  }
}
\details{
  If \code{what="tree"} (the default), the partition tree will be plotted
  using \code{\link[rpart]{plot.rpart}}.

  If \code{what="spatial"}, the predicted intensity
  will be computed using \code{\link{predict.rppm}}, and
  this intensity will be plotted as an image using \code{\link{plot.im}}.
}
\value{
  If \code{what="tree"}, a list containing \code{x} and \code{y}
  coordinates of the plotted nodes of the tree.
  If \code{what="spatial"}, the return value of \code{\link{plot.im}}.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}
}
\seealso{
  \code{\link{rppm}}
}
\examples{
    # Murchison gold data
    mur <- solapply(murchison, rescale, s=1000, unitname="km")
    mur$dfault <- distfun(mur$faults)
    # 
    fit <- rppm(gold ~ dfault + greenstone, data=mur)
    plot(fit)
    plot(fit, what="spatial")
}
\keyword{spatial}
\keyword{hplot}
\keyword{models}
