\name{methods.zclustermodel}
\alias{methods.zclustermodel} % DoNotExport
\alias{pcfmodel.zclustermodel}
\alias{predict.zclustermodel}
\alias{print.zclustermodel}
\title{
  Methods for Cluster Models
}
\description{
  Methods for the experimental class of cluster models.
}
\usage{
 \method{pcfmodel}{zclustermodel}(model, \dots)

 \method{predict}{zclustermodel}(object, \dots,
                  locations, type = "intensity", ngrid = NULL)

 \method{print}{zclustermodel}(x, \dots)
}
\arguments{
  \item{model,object,x}{
    Object of class \code{"zclustermodel"}.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
  \item{locations}{
    Locations where prediction should be performed.
    A window or a point pattern.
  }
  \item{type}{
    Currently must equal \code{"intensity"}.
  }
  \item{ngrid}{
    Pixel grid dimensions for prediction, if \code{locations} is
    a rectangle or polygon.
  }
}
\details{
  Experimental.
}
\value{
  Same as for other methods.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
}
\seealso{
  \code{\link{zclustermodel}}
}
\examples{
  m <- zclustermodel("Thomas", kappa=10, mu=5, scale=0.1)
  m2 <- zclustermodel("VarGamma", kappa=10, mu=10, scale=0.1, nu=0.7)
  m
  m2
  g <- pcfmodel(m)
  g(0.2)
  g2 <- pcfmodel(m2)
  g2(1)
  Z <- predict(m, locations=square(2))
  Z2 <- predict(m2, locations=square(1))
  varcount(m, square(1))
  varcount(m2, square(1))
}
\keyword{spatial}
\keyword{models}
