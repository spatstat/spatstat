\name{G3est}
\Rdversion{1.1}
\alias{G3est}
\title{
  Nearest Neighbour Distance Distribution Function
  of a Three-Dimensional Point Pattern
}
\description{
  Estimates the nearest-neighbour distance distribution function
  \eqn{G_3(r)}{G3(r)} from a three-dimensional point pattern.
}
\usage{
G3est(X, ..., rmax = NULL, nrval = 128, correction = c("rs", "km", "Hanisch"))
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
    \eqn{G_3(r)}{G3(r)} will be estimated. 
  }
  \item{nrval}{
    Optional. Number of values of \eqn{r} for which
    \eqn{G_3(r)}{G3(r)} will be estimated. A large value of \code{nrval}
    is required to avoid discretisation effects.
  }
  \item{correction}{
    Optional. Character vector specifying the edge correction(s)
    to be applied. See Details.
  }
}
\details{
  For a stationary point process \eqn{\Phi}{Phi} in three-dimensional
  space, the nearest-neighbour function
  is
  \deqn{
    G_3(r) = P(d^\ast(x,\Phi) \le r \mid x \in \Phi)
  }{
    G3(r) = P(d*(x,Phi) <= r | x in Phi)
  }
  the cumulative distribution function of the distance 
  \eqn{d^\ast(x,\Phi)}{d*(x,Phi)} from a typical point \eqn{x}
  in  \eqn{\Phi}{Phi} to its nearest neighbour, i.e.
  to the nearest \emph{other} point of \eqn{\Phi}{Phi}.
  
  The three-dimensional point pattern \code{X} is assumed to be a
  partial realisation of a stationary point process \eqn{\Phi}{Phi}.
  The nearest neighbour function of \eqn{\Phi}{Phi} can then be estimated using
  techniques described in the References. For each data point, the
  distance to the nearest neighbour is computed.
  The empirical cumulative distribution
  function of these values, with appropriate edge corrections, is the
  estimate of \eqn{G_3(r)}{G3(r)}.

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
    \item{\code{"Hanisch"}:}{
      the three-dimensional generalisation of
      the Hanisch estimator (Hanisch, 1984).
    }
  }
  Alternatively \code{correction="all"} selects all options.
}
\value{
  A function value table (object of class \code{"fv"}) that can be
  plotted, printed or coerced to a data frame containing the function values.
}
\references{
  Baddeley, A.J, Moyeed, R.A., Howard, C.V. and Boyde, A. (1993)
  Analysis of a three-dimensional point pattern with replication.
  \emph{Applied Statistics} \bold{42}, 641--668.

  Baddeley, A.J. and Gill, R.D. (1997)
  Kaplan-Meier estimators of interpoint distance
  distributions for spatial point processes.
  \emph{Annals of Statistics} \bold{25}, 263--292.

  Hanisch, K.-H. (1984) 
  Some remarks on estimators of the distribution function
  of nearest neighbour distance in stationary spatial point patterns.
  \emph{Mathematische Operationsforschung und Statistik, series Statistics}
  \bold{15}, 409--412.
}
\author{
  \adrian
  and Rana Moyeed.
}
\section{Warnings}{
  A large value of \code{nrval} is required in order to avoid
  discretisation effects (due to the use of histograms in the
  calculation).
}
\seealso{
  \code{\link{pp3}} to create a three-dimensional point
  pattern (object of class \code{"pp3"}).
  
  \code{\link{F3est}},
  \code{\link{K3est}},
  \code{\link{pcf3est}} for other summary functions of
  a three-dimensional point pattern.

  \code{\link{Gest}} to estimate the empty space function of
  point patterns in two dimensions.
}
\examples{
  X <- rpoispp3(42)
  Z <- G3est(X)
  if(interactive()) plot(Z)
}
\keyword{spatial}
\keyword{nonparametric}
