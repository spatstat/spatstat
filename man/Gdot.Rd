\name{Gdot}
\alias{Gdot}
\title{
  Multitype Nearest Neighbour Distance Function (i-to-any)
}
\description{
  For a multitype point pattern, 
  estimate the distribution of the distance
  from a point of type \eqn{i}
  to the nearest other point of any type.
}
\usage{
Gdot(X, i, r=NULL, breaks=NULL, \dots, correction=c("km", "rs", "han"))
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the 
    distance distribution function
    \eqn{G_{i\bullet}(r)}{Gi.(r)} will be computed.
    It must be a multitype point pattern (a marked point pattern
    whose marks are a factor). See under Details.
  }
  \item{i}{The type (mark value)
    of the points in \code{X} from which distances are measured.
    A character string (or something that will be converted to a
    character string).
    Defaults to the first level of \code{marks(X)}.
  }
  \item{r}{Optional. Numeric vector. The values of the argument \eqn{r}
    at which the distribution function
    \eqn{G_{i\bullet}(r)}{Gi.(r)} should be evaluated.
    There is a sensible default.
    First-time users are strongly advised not to specify this argument.
    See below for important conditions on \eqn{r}.
  }
  \item{breaks}{
    This argument is for internal use only.
  }
  \item{\dots}{Ignored.}
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
    at which the function \eqn{G_{i\bullet}(r)}{Gi.(r)} has been  estimated
  }
  \item{rs}{the ``reduced sample'' or ``border correction''
    estimator of \eqn{G_{i\bullet}(r)}{Gi.(r)}
  }
  \item{han}{the Hanisch-style estimator of \eqn{G_{i\bullet}(r)}{Gi.(r)}
  }
  \item{km}{the spatial Kaplan-Meier estimator of \eqn{G_{i\bullet}(r)}{Gi.(r)}
  }
  \item{hazard}{the hazard rate \eqn{\lambda(r)}{lambda(r)}
    of \eqn{G_{i\bullet}(r)}{Gi.(r)} by the spatial Kaplan-Meier method
  }
  \item{raw}{the uncorrected estimate of \eqn{G_{i\bullet}(r)}{Gi.(r)},
  i.e. the empirical distribution of the distances from 
  each point of type \eqn{i} to the nearest other point of any type.
  }
  \item{theo}{the theoretical value of \eqn{G_{i\bullet}(r)}{Gi.(r)}
    for a marked Poisson process with the same estimated intensity
    (see below).
  }
}
\details{
  This function \code{Gdot} and its companions
  \code{\link{Gcross}} and \code{\link{Gmulti}}
  are generalisations of the function \code{\link{Gest}}
  to multitype point patterns.

  A multitype point pattern is a spatial pattern of
  points classified into a finite number of possible
  ``colours'' or ``types''. In the \pkg{spatstat} package,
  a multitype pattern is represented as a single 
  point pattern object in which the points carry marks,
  and the mark value attached to each point
  determines the type of that point.
  
  The argument \code{X} must be a point pattern (object of class
  \code{"ppp"}) or any data that are acceptable to \code{\link{as.ppp}}.
  It must be a marked point pattern, and the mark vector
  \code{X$marks} must be a factor.
  The argument will be interpreted as a
  level of the factor \code{X$marks}. (Warning: this means that
  an integer value \code{i=3} will be interpreted as the number 3,
  \bold{not} the 3rd smallest level.)
  
  The ``dot-type'' (type \eqn{i} to any type)
  nearest neighbour distance distribution function 
  of a multitype point process 
  is the cumulative distribution function \eqn{G_{i\bullet}(r)}{Gi.(r)}
  of the distance from a typical random point of the process with type \eqn{i}
  the nearest other point of the process, regardless of type. 

  An estimate of \eqn{G_{i\bullet}(r)}{Gi.(r)}
  is a useful summary statistic in exploratory data analysis
  of a multitype point pattern.
  If the type \eqn{i} points
  were independent of all other points,
  then \eqn{G_{i\bullet}(r)}{Gi.(r)} would equal \eqn{G_{ii}(r)}{Gii(r)},
  the nearest neighbour distance distribution function of the type
  \eqn{i} points alone.
  For a multitype Poisson point process with total intensity
  \eqn{\lambda}{lambda}, we have
  \deqn{G_{i\bullet}(r) = 1 - e^{ - \lambda \pi r^2} }{%
    Gi.(r) = 1 - exp( - lambda * pi * r^2)}
  Deviations between the empirical and theoretical
  \eqn{G_{i\bullet}}{Gi.} curves
  may suggest dependence of the type \eqn{i} points on the other points.

  This algorithm estimates the distribution function
  \eqn{G_{i\bullet}(r)}{Gi.(r)} 
  from the point pattern \code{X}. It assumes that \code{X} can be treated
  as a realisation of a stationary (spatially homogeneous) 
  random spatial point process in the plane, observed through
  a bounded window.
  The window (which is specified in \code{X} as \code{Window(X)})
  may have arbitrary shape.
  Biases due to edge effects are
  treated in the same manner as in \code{\link{Gest}}.

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{G_{i\bullet}(r)}{Gi.(r)} should be evaluated. 
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
  \eqn{\lambda(r)}{lambda(r)}, of \eqn{G_{i\bullet}(r)}{Gi.(r)}. 
  This estimate should be used with caution as
  \eqn{G_{i\bullet}(r)}{Gi.(r)}
  is not necessarily differentiable.

  The naive empirical distribution of distances from each point of
  the pattern \code{X} to the nearest other point of the pattern, 
  is a biased estimate of \eqn{G_{i\bullet}}{Gi.}.
  However this is also returned by the algorithm, as it is sometimes 
  useful in other contexts. Care should be taken not to use the uncorrected
  empirical \eqn{G_{i\bullet}}{Gi.} as if it were an unbiased estimator of
  \eqn{G_{i\bullet}}{Gi.}.
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
  The argument \code{i} is interpreted as
  a level of the factor \code{X$marks}. It is converted to a character
  string if it is not already a character string.
  The value \code{i=1} does \bold{not}
  refer to the first level of the factor.

  The function \eqn{G_{i\bullet}}{Gi.} does not necessarily have a density. 

  The reduced sample estimator of \eqn{G_{i\bullet}}{Gi.}
  is pointwise approximately 
  unbiased, but need not be a valid distribution function; it may 
  not be a nondecreasing function of \eqn{r}. Its range is always 
  within \eqn{[0,1]}.

  The spatial Kaplan-Meier estimator of \eqn{G_{i\bullet}}{Gi.}
  is always nondecreasing
  but its maximum value may be less than \eqn{1}.
}
\seealso{
 \code{\link{Gcross}},
 \code{\link{Gest}},
 \code{\link{Gmulti}}
}
\examples{
    # amacrine cells data
    G0. <- Gdot(amacrine, "off") 
    plot(G0.)

    # synthetic example    
    pp <- runifpoispp(30)
    pp <- pp \%mark\% factor(sample(0:1, npoints(pp), replace=TRUE))
    G <- Gdot(pp, "0")
    G <- Gdot(pp, 0) # equivalent
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{nonparametric}




