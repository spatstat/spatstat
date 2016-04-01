\name{Jmulti}
\alias{Jmulti}
\title{
  Marked J Function
}
\description{
  For a marked point pattern, 
  estimate the multitype \eqn{J} function
  summarising dependence between the
  points in subset \eqn{I}
  and those in subset \eqn{J}.
}
\usage{
  Jmulti(X, I, J, eps=NULL, r=NULL, breaks=NULL, \dots, disjoint=NULL,
         correction=NULL)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the multitype distance distribution function
    \eqn{J_{IJ}(r)}{J[IJ](r)} will be computed.
    It must be a marked point pattern.
    See under Details.
  }
  \item{I}{Subset of points of \code{X} from which distances are
    measured. See Details.
  }
  \item{J}{Subset of points in \code{X} to which distances are measured.
    See Details.
  }
  \item{eps}{A positive number.
    The pixel resolution of the discrete approximation to Euclidean
    distance (see \code{\link{Jest}}). There is a sensible default.
  }
  \item{r}{numeric vector. The values of the argument \eqn{r}
    at which the distribution function
    \eqn{J_{IJ}(r)}{J[IJ](r)} should be evaluated.
    There is a sensible default.
    First-time users are strongly advised not to specify this argument.
    See below for important conditions on \eqn{r}.
  }
  \item{breaks}{
	This argument is for internal use only.
  }
  \item{\dots}{Ignored.}
  \item{disjoint}{Optional flag indicating whether
    the subsets \code{I} and \code{J} are disjoint.
    If missing, this value will be computed by inspecting the
    vectors \code{I} and \code{J}.
  }
  \item{correction}{
    Optional. Character string specifying the edge correction(s)
    to be used. Options are \code{"none"}, \code{"rs"}, \code{"km"},
    \code{"Hanisch"} and \code{"best"}.
    Alternatively \code{correction="all"} selects all options.
  }
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).

  Essentially a data frame containing six numeric columns 
  \item{r}{the values of the argument \eqn{r} 
    at which the function \eqn{J_{IJ}(r)}{J[IJ](r)} has been  estimated
  }
  \item{rs}{the ``reduced sample'' or ``border correction''
    estimator of \eqn{J_{IJ}(r)}{J[IJ](r)}
  }
  \item{km}{the spatial Kaplan-Meier estimator of \eqn{J_{IJ}(r)}{J[IJ](r)}
  }
  \item{han}{the Hanisch-style estimator of \eqn{J_{IJ}(r)}{J[IJ](r)}
  }
  \item{un}{the uncorrected estimate of \eqn{J_{IJ}(r)}{J[IJ](r)},
    formed by taking the ratio of uncorrected empirical estimators
    of \eqn{1 - G_{IJ}(r)}{1 - G[IJ](r)}
    and \eqn{1 - F_{J}(r)}{1 - F[J](r)}, see
    \code{\link{Gdot}} and \code{\link{Fest}}.
  }
  \item{theo}{the theoretical value of \eqn{J_{IJ}(r)}{J[IJ](r)}
    for a marked Poisson process with the same estimated intensity,
    namely 1.
  }
}
\details{
  The function \code{Jmulti}
  generalises \code{\link{Jest}} (for unmarked point
  patterns) and \code{\link{Jdot}} and \code{\link{Jcross}} (for
  multitype point patterns) to arbitrary marked point patterns.

  Suppose \eqn{X_I}{X[I]}, \eqn{X_J}{X[J]} are subsets, possibly
  overlapping, of a marked point process. Define
  \deqn{J_{IJ}(r) = \frac{1 - G_{IJ}(r)}{1 - F_J(r)}}{
    J[IJ](r) = (1 - G[IJ](r))/(1 - F[J](r))}
  where \eqn{F_J(r)}{F[J](r)} is the cumulative distribution function of
  the distance from a fixed location to the nearest point
  of \eqn{X_J}{X[J]}, and \eqn{G_{IJ}(r)}{GJ(r)}
  is the distribution function of the distance
  from a typical point of  \eqn{X_I}{X[I]} to the nearest distinct point of
  \eqn{X_J}{X[J]}. 

  The argument \code{X} must be a point pattern (object of class
  \code{"ppp"}) or any data that are acceptable to \code{\link{as.ppp}}.

  The arguments \code{I} and \code{J} specify two subsets of the
  point pattern. They may be any type of subset indices, for example,
  logical vectors of length equal to \code{npoints(X)},
  or integer vectors with entries in the range 1 to
  \code{npoints(X)}, or negative integer vectors.

  Alternatively, \code{I} and \code{J} may be \bold{functions}
  that will be applied to the point pattern \code{X} to obtain
  index vectors. If \code{I} is a function, then evaluating
  \code{I(X)} should yield a valid subset index. This option
  is useful when generating simulation envelopes using
  \code{\link{envelope}}.

  It is assumed that \code{X} can be treated
  as a realisation of a stationary (spatially homogeneous) 
  random spatial point process in the plane, observed through
  a bounded window.
  The window (which is specified in \code{X} as \code{Window(X)})
  may have arbitrary shape.
  Biases due to edge effects are
  treated in the same manner as in \code{\link{Jest}}.

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{J_{IJ}(r)}{J[IJ](r)} should be evaluated. 
  It is also used to determine the breakpoints
  (in the sense of \code{\link{hist}})
  for the computation of histograms of distances. The reduced-sample and
  Kaplan-Meier estimators are computed from histogram counts. 
  In the case of the Kaplan-Meier estimator this introduces a discretisation
  error which is controlled by the fineness of the breakpoints.

  First-time users would be strongly advised not to specify \code{r}.
  However, if it is specified, \code{r} must satisfy \code{r[1] = 0}, 
  and \code{max(r)} must be larger than the radius of the largest disc 
  contained in the window. Furthermore, the successive entries of \code{r}
  must be finely spaced.
}
\references{
  Van Lieshout, M.N.M. and Baddeley, A.J. (1999)
  Indices of dependence between types in multivariate point patterns.
  \emph{Scandinavian Journal of Statistics} \bold{26}, 511--532.

}
\seealso{
 \code{\link{Jcross}},
 \code{\link{Jdot}},
 \code{\link{Jest}}
}
\examples{
    trees <- longleaf
     # Longleaf Pine data: marks represent diameter
    \testonly{
        trees <- trees[seq(1,npoints(trees), by=50)]
    }
    Jm <- Jmulti(trees, marks(trees) <= 15, marks(trees) >= 25)
    plot(Jm)
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}.
}
\keyword{spatial}
\keyword{nonparametric}

