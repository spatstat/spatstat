\name{F3est}
\Rdversion{1.1}
\alias{F3est}
\title{
  Empty Space Function of a Three-Dimensional Point Pattern
}
\description{
  Estimates the empty space function \eqn{F_3(r)}{F3(r)} from 
  a three-dimensional point pattern.
}
\usage{
F3est(X, ..., rmax = NULL, nrval = 128, vside = NULL,
              correction = c("rs", "km", "cs"),
              sphere = c("fudge", "ideal", "digital"))
}
\arguments{
  \item{X}{
    Three-dimensional point pattern (object of class \code{"pp3"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{rmax}{
    Optional. Maximum value of argument \eqn{r} for which
    \eqn{F_3(r)}{F3(r)} will be estimated. 
  }
  \item{nrval}{
    Optional. Number of values of \eqn{r} for which
    \eqn{F_3(r)}{F3(r)} will be estimated. A large value of \code{nrval}
    is required to avoid discretisation effects.
  }
  \item{vside}{
    Optional. 
    Side length of the voxels in the discrete approximation.
  }
  \item{correction}{
    Optional. Character vector specifying the edge correction(s)
    to be applied. See Details.
  }
  \item{sphere}{
    Optional. Character string specifying how to calculate the
    theoretical value of \eqn{F_3(r)}{F3(r)} for a Poisson
    process. See Details.
  }
}
\details{
  For a stationary point process \eqn{\Phi}{Phi} in three-dimensional
  space, the empty space function is
  \deqn{
    F_3(r) = P(d(0,\Phi) \le r)
  }{
    F3(r) = P(d(0,Phi) <= r)
  }
  where \eqn{d(0,\Phi)}{d(0,Phi)} denotes the distance from a fixed
  origin \eqn{0} to the nearest point of \eqn{\Phi}{Phi}.
  
  The three-dimensional point pattern \code{X} is assumed to be a
  partial realisation of a stationary point process \eqn{\Phi}{Phi}.
  The empty space function of \eqn{\Phi}{Phi} can then be estimated using
  techniques described in the References.

  The box containing the point
  pattern is discretised into cubic voxels of side length \code{vside}.
  The distance function \eqn{d(u,\Phi)}{d(u,Phi)} is computed for
  every voxel centre point
  \eqn{u} using a three-dimensional version of the distance transform
  algorithm (Borgefors, 1986). The empirical cumulative distribution
  function of these values, with appropriate edge corrections, is the
  estimate of \eqn{F_3(r)}{F3(r)}.

  The available edge corrections are:
  \describe{
    \item{\code{"rs"}:}{
      the reduced sample (aka minus sampling, border correction)
      estimator (Baddeley et al, 1993)
    }
    \item{\code{"km"}:}{
      the three-dimensional version of the
      Kaplan-Meier estimator (Baddeley and Gill, 1997)
    }
    \item{\code{"cs"}:}{
      the three-dimensional generalisation of
      the Chiu-Stoyan or Hanisch estimator (Chiu and Stoyan, 1998).
    }
  }

  The result includes a column \code{theo} giving the 
  theoretical value of \eqn{F_3(r)}{F3(r)} for
  a uniform Poisson process (Complete Spatial Randomness).
  This value depends on the volume of the sphere of radius \code{r}
  measured in the discretised distance metric.
  The argument \code{sphere} determines how this will be calculated.
  \itemize{
    \item
    If \code{sphere="ideal"} the calculation will use the
    volume of an ideal sphere of radius \eqn{r} namely
    \eqn{(4/3) \pi r^3}{(4/3) * pi * r^3}. This is not recommended
    because the theoretical values of \eqn{F_3(r)}{F3(r)}
    are inaccurate.
    \item
    If \code{sphere="fudge"} then the volume of the ideal sphere will
    be multiplied by 0.78, which gives the approximate volume
    of the sphere in the discretised distance metric.
    \item
    If \code{sphere="digital"} then the volume of the sphere in the
    discretised distance metric is computed exactly using another
    distance transform. This takes longer to compute, but is exact.
  }
}
\value{
  A function value table (object of class \code{"fv"}) that can be
  plotted, printed or coerced to a data frame containing the function values.
}
\references{
  Baddeley, A.J, Moyeed, R.A., Howard, C.V. and Boyde, A.
  Analysis of a three-dimensional point pattern with replication.
  \emph{Applied Statistics} \bold{42} (1993) 641--668.

  Baddeley, A.J. and Gill, R.D. (1997)
  Kaplan-Meier estimators of interpoint distance
  distributions for spatial point processes.
  \emph{Annals of Statistics} \bold{25}, 263--292.

  Borgefors, G. (1986)
  Distance transformations in digital images.
  \emph{Computer Vision, Graphics and Image Processing}
  \bold{34}, 344--371.

  Chiu, S.N. and Stoyan, D. (1998)
  Estimators of distance distributions for spatial patterns.
  \emph{Statistica Neerlandica} \bold{52}, 239--246.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rana Moyeed.
}
\section{Warnings}{
  A small value of \code{vside} and a large value of \code{nrval}
  are required for reasonable accuracy. 

  The default value of \code{vside} ensures that the total number of
  voxels is \code{2^22} or about 4 million.
  To change the default number of voxels, see
  \code{\link{spatstat.options}("nvoxel")}.
}
\seealso{
  \code{\link{G3est}},
  \code{\link{K3est}},
  \code{\link{pcf3est}}.
}
\examples{
  \testonly{op <- spatstat.options(nvoxel=2^18)}
  X <- rpoispp3(42)
  Z <- F3est(X)
  if(interactive()) plot(Z)
  \testonly{spatstat.options(op)}
}
\keyword{spatial}
\keyword{nonparametric}
