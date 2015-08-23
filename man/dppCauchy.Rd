\name{dppCauchy}
\alias{dppCauchy}
\title{Generalized Cauchy Determinantal Point Process Model}
\description{
  Function generating an instance of the (generalized) Cauchy
  determinantal point process model.
}
\usage{dppCauchy(\dots)}
\arguments{
  \item{\dots}{arguments of the form \code{tag=value} specifying the
    parameters. See Details.} 
}
\details{
  The possible parameters are:
  \itemize{
    \item the intensity \code{lambda} as a positive numeric
    \item the scale parameter \code{alpha} as a positive numeric
    \item the shape parameter \code{nu} as a positive numeric
    (artificially required to be less than 20 in the code for numerical
    stability)
    \item the dimension \code{d} as a positive integer
  }
}
\value{An object of class \code{"dppmodel"}.}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\examples{
m <- dppCauchy(lambda=100, alpha=.05, nu=1, d=2)
}
\seealso{
  \code{\link{dppBessel}},
  \code{\link{dppGauss}},
  \code{\link{dppMatern}},
  \code{\link{dppPowerExp}}
}

