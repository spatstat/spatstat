\name{Gcross}
\alias{Gcross}
\title{
  Multitype Nearest Neighbour Distance Function (i-to-j)
}
\description{
  For a multitype point pattern, 
  estimate the distribution of the distance
  from a point of type \eqn{i}
  to the nearest point of type \eqn{j}.
}
\usage{
Gcross(X, i, j, r=NULL, breaks=NULL, \dots, correction=c("rs", "km", "han"))
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the cross type distance distribution function
    \eqn{G_{ij}(r)}{Gij(r)} will be computed.
    It must be a multitype point pattern (a marked point pattern
    whose marks are a factor). See under Details.
  }
  \item{i}{The type (mark value)
    of the points in \code{X} from which distances are measured.
    A character string (or something that will be converted to a
    character string).
    Defaults to the first level of \code{marks(X)}.
  }
  \item{j}{The type (mark value)
    of the points in \code{X} to which distances are measured.
    A character string (or something that will be
    converted to a character string).
    Defaults to the second level of \code{marks(X)}.
  }
  \item{r}{Optional. Numeric vector. The values of the argument \eqn{r}
    at which the distribution function
    \eqn{G_{ij}(r)}{Gij(r)} should be evaluated.
    There is a sensible default.
    First-time users are strongly advised not to specify this argument.
    See below for important conditions on \eqn{r}.
  }
  \item{breaks}{An alternative to the argument \code{r}.
    Not normally invoked by the user. See the \bold{Details} section.
  }
  \item{\dots}{
  	Ignored.
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
    at which the function \eqn{G_{ij}(r)}{Gij(r)} has been  estimated
  }
  \item{rs}{the ``reduced sample'' or ``border correction''
    estimator of \eqn{G_{ij}(r)}{Gij(r)}
  }
  \item{han}{the Hanisch-style estimator of \eqn{G_{ij}(r)}{Gij(r)}
  }
  \item{km}{the spatial Kaplan-Meier estimator of \eqn{G_{ij}(r)}{Gij(r)}
  }
  \item{hazard}{the hazard rate \eqn{\lambda(r)}{lambda(r)}
    of \eqn{G_{ij}(r)}{Gij(r)} by the spatial Kaplan-Meier method
  }
  \item{raw}{the uncorrected estimate of \eqn{G_{ij}(r)}{Gij(r)},
  i.e. the empirical distribution of the distances from 
  each point of type \eqn{i} to the nearest point of type \eqn{j}
  }
  \item{theo}{the theoretical value of \eqn{G_{ij}(r)}{Gij(r)}
    for a marked Poisson process with the same estimated intensity
    (see below).
  }
}
\details{
  This function \code{Gcross} and its companions
  \code{\link{Gdot}} and \code{\link{Gmulti}}
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
  The arguments \code{i} and \code{j} will be interpreted as
  levels of the factor \code{X$marks}. (Warning: this means that
  an integer value \code{i=3} will be interpreted as
  the number 3, \bold{not} the 3rd smallest level). 
  
  The ``cross-type'' (type \eqn{i} to type \eqn{j})
  nearest neighbour distance distribution function 
  of a multitype point process 
  is the cumulative distribution function \eqn{G_{ij}(r)}{Gij(r)}
  of the distance from a typical random point of the process with type \eqn{i}
  the nearest point of type \eqn{j}. 

  An estimate of \eqn{G_{ij}(r)}{Gij(r)}
  is a useful summary statistic in exploratory data analysis
  of a multitype point pattern.
  If the process of type \eqn{i} points
  were independent of the process of type \eqn{j} points,
  then \eqn{G_{ij}(r)}{Gij(r)} would equal \eqn{F_j(r)}{Fj(r)},
  the empty space function of the type \eqn{j} points.
  For a multitype Poisson point process where the type \eqn{i} points
  have intensity \eqn{\lambda_i}{lambda[i]}, we have
  \deqn{G_{ij}(r) = 1 - e^{ - \lambda_j \pi r^2} }{%
    Gij(r) = 1 - exp( - lambda[j] * pi * r^2)}
  Deviations between the empirical and theoretical \eqn{G_{ij}}{Gij} curves
  may suggest dependence between the points of types \eqn{i} and \eqn{j}.

  This algorithm estimates the distribution function \eqn{G_{ij}(r)}{Gij(r)} 
  from the point pattern \code{X}. It assumes that \code{X} can be treated
  as a realisation of a stationary (spatially homogeneous) 
  random spatial point process in the plane, observed through
  a bounded window.
  The window (which is specified in \code{X} as \code{X$window})
  may have arbitrary shape.
  Biases due to edge effects are
  treated in the same manner as in \code{\link{Gest}}.

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{G_{ij}(r)}{Gij(r)} should be evaluated. 
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
  \eqn{\lambda(r)}{lambda(r)}, of \eqn{G_{ij}(r)}{Gij(r)}. 
  This estimate should be used with caution as \eqn{G_{ij}(r)}{Gij(r)}
  is not necessarily differentiable.

  The naive empirical distribution of distances from each point of
  the pattern \code{X} to the nearest other point of the pattern, 
  is a biased estimate of \eqn{G_{ij}}{Gij}.
  However this is also returned by the algorithm, as it is sometimes 
  useful in other contexts. Care should be taken not to use the uncorrected
  empirical \eqn{G_{ij}}{Gij} as if it were an unbiased estimator of
  \eqn{G_{ij}}{Gij}.
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
  The arguments \code{i} and \code{j} are always interpreted as
  levels of the factor \code{X$marks}. They are converted to character
  strings if they are not already character strings.
  The value \code{i=1} does \bold{not}
  refer to the first level of the factor.

  The function \eqn{G_{ij}}{Gij} does not necessarily have a density. 

  The reduced sample estimator of \eqn{G_{ij}}{Gij} is pointwise approximately 
  unbiased, but need not be a valid distribution function; it may 
  not be a nondecreasing function of \eqn{r}. Its range is always 
  within \eqn{[0,1]}.

  The spatial Kaplan-Meier estimator of \eqn{G_{ij}}{Gij}
  is always nondecreasing
  but its maximum value may be less than \eqn{1}.
}
\seealso{
 \code{\link{Gdot}},
 \code{\link{Gest}},
 \code{\link{Gmulti}}
}
\examples{
    # amacrine cells data
    G01 <- Gcross(amacrine)

    # equivalent to:
    \dontrun{
    G01 <- Gcross(amacrine, "off", "on")
    }

    plot(G01)

    # empty space function of `on' points
    \dontrun{
       F1 <- Fest(split(amacrine)$on, r = G01$r)
       lines(F1$r, F1$km, lty=3)
    }

    # synthetic example    
    pp <- runifpoispp(30)
    pp <- pp \%mark\% factor(sample(0:1, pp$n, replace=TRUE))
    G <- Gcross(pp, "0", "1")   # note: "0" not 0
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
