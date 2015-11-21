\name{zclustermodel}
\alias{zclustermodel}
\title{
  Cluster Point Process Model
}
\description{
  Experimental code. Creates an object representing a cluster point
  process model. Typically used for theoretical calculations about
  such a model.
}
\usage{
zclustermodel(name = "Thomas", \dots, mu, kappa, scale)
}
\arguments{
  \item{name}{
    Name of the cluster process. One of
    \code{"Thomas"},
    \code{"MatClust"},
    \code{"VarGamma"} or
    \code{"Cauchy"}.
  }
  \item{\dots}{
    Other arguments needed for the model.
  }
  \item{mu}{
    Mean cluster size. A single number, or a pixel image.
  }
  \item{kappa}{
    Parent intensity. A single number.
  }
  \item{scale}{
    Cluster scale parameter of the model.
  }
}
\details{
  Experimental.
}
\value{
  Object of the experimental class \code{"zclustermodel"}.
}
\author{
  \adrian
}
\seealso{
  \code{\link{methods.zclustermodel}}
}
\examples{
  m <- zclustermodel("Thomas", kappa=10, mu=5, scale=0.1)
}
\keyword{spatial}
\keyword{models}

