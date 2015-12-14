\name{clusterradius}
\alias{clusterradius}
\alias{clusterradius.character}
\alias{clusterradius.kppm}
\title{
  Compute or Extract Effective Range of Cluster Kernel
}
\description{
  Given a cluster point process model, this command
  returns a value beyond which the the probability density of the
  cluster offspring is neglible.
}
\usage{
clusterradius(model, \dots)

\method{clusterradius}{kppm}(model, \dots, thresh = NULL, precision = FALSE)

\method{clusterradius}{character}(model, \dots, thresh = NULL, precision = FALSE)
}
\arguments{
  \item{model}{
    Cluster model. Either a fitted cluster or Cox model
    (object of class \code{"kppm"}), or a character string
    specifying the type of cluster model.
  }
  \item{\dots}{
    Parameter values for the model,
    when \code{model} is a character string.
  }
  \item{thresh}{
    Numerical threshold relative to the cluster kernel value at the
    origin (parent location) determining when the cluster kernel
    will be considered neglible. A sensible default is provided.
  }
  \item{precision}{
    Logical. If \code{precision=TRUE} the precision of the calculated
    range is returned as an attribute to the range. See details.
  }
}
\details{
  Given a cluster model this function by default returns the effective
  range of the model with the given parameters as used in spatstat. For
  the Matern cluster model (see e.g. \code{\link{rMatClust}}) this is
  simply the finite radius of the offsring density given by the paramter
  \code{scale} irrespective of other options given to this function. The
  remaining models in spatstat have infinite theoretical range, and an
  effective finite value is given as follows: For the Thomas model (see
  e.g. \code{\link{rThomas}} the default is \code{4*scale} where scale
  is the scale or standard deviation parameter of the model. If
  \code{thresh} is given the value is instead found as described for the
  other models below.

  For the Cauchy model (see e.g. \code{\link{rCauchy}}) and the Variance
  Gamma (Bessel) model (see e.g. \code{\link{rVarGamma}}) the value of
  \code{thresh} defaults to 0.001, and then this is used to compute the
  range numerically as follows. If \eqn{k(x,y)=k_0(r)} with \eqn{r=\sqrt(x^2+y^2)}{r=sqrt(x^2+y^2)}
  denotes the isotropic cluster kernel then \eqn{f(r) = 2 \pi r
  k_0(r)}{f(r) = 2 pi r k_0(r)} is the
  density function of the offspring distance from the parent. The range
  is determined as the value of \eqn{r} where \eqn{f(r)} falls below
  \code{thresh} times \eqn{k_0(r)}.
  
  If \code{precision=TRUE} the precision related to the chosen range is
  returned as an attribute. Here the precision is defined as the polar
  integral of the kernel from distance 0 to the calculated
  range. Ideally this should be close to the value 1 which would be
  obtained for the true theretical infinite range.

}
\value{
  A positive numeric.

  Additionally, the precision related to this range value is returned as
  an attribute \code{"prec"}, if \code{precision=TRUE}.
}
\author{
Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}


Rolf Turner \email{r.turner@auckland.ac.nz}

and Ege Rubak \email{rubak@math.aau.dk}

}
\seealso{
  \code{\link{clusterkernel}}, \code{\link{kppm}},
  \code{\link{rMatClust}}, \code{\link{rThomas}}, \code{\link{rCauchy}},
  \code{\link{rVarGamma}}, \code{\link{rNeymanScott}}.
}
\examples{
  fit <- kppm(redwood ~ x, "MatClust")
  clusterradius(fit)

  clusterradius("Thomas", scale = .1)
  clusterradius("Thomas", scale = .1, thresh = 0.001)
  clusterradius("VarGamma", scale = .1, nu = 2, precision = TRUE)
}
\keyword{spatial}
