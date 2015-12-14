\name{dppBessel}
\alias{dppBessel}
\title{Bessel Type Determinantal Point Process Model}
\description{
  Function generating an instance of the Bessel-type
  determinantal point process model.
}
\usage{dppBessel(\dots)}
\arguments{
  \item{\dots}{arguments of the form \code{tag=value}
    specifying the model parameters. See Details.
  }
}
\details{
  The possible parameters are:
  \itemize{
    \item the intensity \code{lambda} as a positive numeric
    \item the scale parameter \code{alpha} as a positive numeric
    \item the shape parameter \code{sigma} as a non-negative numeric
    \item the dimension \code{d} as a positive integer
  }
}
\value{An object of class \code{"detpointprocfamily"}.}
\author{
  Frederic Lavancier and Christophe Biscio.
  Modified by Ege Rubak \email{rubak@math.aau.dk} ,
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  and 
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\examples{
m <- dppBessel(lambda=100, alpha=.05, sigma=0, d=2)
}
\seealso{
  \code{\link{dppCauchy}},
  \code{\link{dppGauss}},
  \code{\link{dppMatern}},
  \code{\link{dppPowerExp}}
}

