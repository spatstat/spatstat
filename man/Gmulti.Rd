\name{Gmulti}
\alias{Gmulti}
\title{
  Marked Nearest Neighbour Distance Function
}
\description{
  For a marked point pattern, 
  estimate the distribution of the distance
  from a typical point in subset \code{I}
  to the nearest point of subset \eqn{J}.
}
\usage{
Gmulti(X, I, J, r=NULL, breaks=NULL, \dots,
        disjoint=NULL, correction=c("rs", "km", "han"))
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the multitype distance distribution function
    \eqn{G_{IJ}(r)}{GIJ(r)} will be computed.
    It must be a marked point pattern.
    See under Details.
  }
  \item{I}{Subset of points of \code{X} from which distances are
    measured. 
  }
  \item{J}{Subset of points in \code{X} to which distances are measured.
  }
  \item{r}{Optional. Numeric vector. The values of the argument \eqn{r}
    at which the distribution function
    \eqn{G_{IJ}(r)}{GIJ(r)} should be evaluated.
    There is a sensible default.
    First-time users are strongly advised not to specify this argument.
    See below for important conditions on \eqn{r}.
  }
  \item{breaks}{An alternative to the argument \code{r}.
    Not normally invoked by the user. See the \bold{Details} section.
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
    \code{"hanisch"} and \code{"best"}.
  }
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).
  
  Essentially a data frame containing six numeric columns 
  \item{r}{the values of the argument \eqn{r} 
    at which the function \eqn{G_{IJ}(r)}{GIJ(r)} has been  estimated
  }
  \item{rs}{the ``reduced sample'' or ``border correction''
    estimator of \eqn{G_{IJ}(r)}{GIJ(r)}
  }
  \item{han}{the Hanisch-style estimator of \eqn{G_{IJ}(r)}{GIJ(r)}
  }
  \item{km}{the spatial Kaplan-Meier estimator of \eqn{G_{IJ}(r)}{GIJ(r)}
  }
  \item{hazard}{the hazard rate \eqn{\lambda(r)}{lambda(r)}
    of \eqn{G_{IJ}(r)}{GIJ(r)} by the spatial Kaplan-Meier method
  }
  \item{raw}{the uncorrected estimate of \eqn{G_{IJ}(r)}{GIJ(r)},
  i.e. the empirical distribution of the distances from 
  each point of type \eqn{i} to the nearest point of type \eqn{j}
  }
  \item{theo}{the theoretical value of \eqn{G_{IJ}(r)}{GIJ(r)}
    for a marked Poisson process with the same estimated intensity
  }
}
\details{
  The function \code{Gmulti}
  generalises \code{\link{Gest}} (for unmarked point
  patterns) and \code{\link{Gdot}} and \code{\link{Gcross}} (for
  multitype point patterns) to arbitrary marked point patterns.

  Suppose \eqn{X_I}{X[I]}, \eqn{X_J}{X[J]} are subsets, possibly
  overlapping, of a marked point process. This function computes an
  estimate of the cumulative
  distribution function \eqn{G_{IJ}(r)}{GIJ(r)} of the distance
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

  This algorithm estimates the distribution function \eqn{G_{IJ}(r)}{GIJ(r)} 
  from the point pattern \code{X}. It assumes that \code{X} can be treated
  as a realisation of a stationary (spatially homogeneous) 
  random spatial point process in the plane, observed through
  a bounded window.
  The window (which is specified in \code{X} as \code{X$window})
  may have arbitrary shape.
  Biases due to edge effects are
  treated in the same manner as in \code{\link{Gest}}.

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{G_{IJ}(r)}{GIJ(r)} should be evaluated. 
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

  The algorithm also returns an estimate of the hazard rate function, 
  \eqn{\lambda(r)}{lambda(r)}, of \eqn{G_{IJ}(r)}{GIJ(r)}. 
  This estimate should be used with caution as \eqn{G_{IJ}(r)}{GIJ(r)}
  is not necessarily differentiable.

  The naive empirical distribution of distances from each point of
  the pattern \code{X} to the nearest other point of the pattern, 
  is a biased estimate of \eqn{G_{IJ}}{GIJ}.
  However this is also returned by the algorithm, as it is sometimes 
  useful in other contexts. Care should be taken not to use the uncorrected
  empirical \eqn{G_{IJ}}{GIJ} as if it were an unbiased estimator of
  \eqn{G_{IJ}}{GIJ}.
}
\references{
  Cressie, N.A.C. \emph{Statistics for spatial data}.
    John Wiley and Sons, 1991.

  Diggle, P.J. \emph{Statistical analysis of spatial point patterns}.
    Academic Press, 1983.

  Diggle, P. J. (1986).
  Displaced amacrine cells in the retina of a
  rabbit : analysis of a bivariate spatial point pattern. 
  \emph{J. Neurosci. Meth.} \bold{18}, 115--125.
 
  Harkness, R.D and Isham, V. (1983)
  A bivariate spatial point pattern of ants' nests.
  \emph{Applied Statistics} \bold{32}, 293--303
 
  Lotwick, H. W. and Silverman, B. W. (1982).
  Methods for analysing spatial processes of several types of points.
  \emph{J. Royal Statist. Soc. Ser. B} \bold{44}, 406--413.

  Ripley, B.D. \emph{Statistical inference for spatial processes}.
  Cambridge University Press, 1988.

  Stoyan, D, Kendall, W.S. and Mecke, J.
  \emph{Stochastic geometry and its applications}.
  2nd edition. Springer Verlag, 1995.

  Van Lieshout, M.N.M. and Baddeley, A.J. (1999)
  Indices of dependence between types in multivariate point patterns.
  \emph{Scandinavian Journal of Statistics} \bold{26}, 511--532.

}
\section{Warnings}{
  The function \eqn{G_{IJ}}{GIJ} does not necessarily have a density. 

  The reduced sample estimator of \eqn{G_{IJ}}{GIJ} is pointwise approximately 
  unbiased, but need not be a valid distribution function; it may 
  not be a nondecreasing function of \eqn{r}. Its range is always 
  within \eqn{[0,1]}.

  The spatial Kaplan-Meier estimator of \eqn{G_{IJ}}{GIJ}
  is always nondecreasing
  but its maximum value may be less than \eqn{1}.
}
\seealso{
 \code{\link{Gcross}},
 \code{\link{Gdot}},
 \code{\link{Gest}}
}
\examples{
    data(longleaf)
     # Longleaf Pine data: marks represent diameter
    \testonly{
      longleaf <- longleaf[seq(1,longleaf$n, by=50), ]
    }
    Gm <- Gmulti(longleaf, longleaf$marks <= 15, longleaf$marks >= 25)
    plot(Gm)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}

