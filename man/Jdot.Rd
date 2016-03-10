\name{Jdot}
\alias{Jdot}
\title{
  Multitype J Function (i-to-any)
}
\description{
  For a multitype point pattern, 
  estimate the multitype \eqn{J} function 
  summarising the interpoint dependence between
  the type \eqn{i} points and the points of any type.
}
\usage{
Jdot(X, i, eps=NULL, r=NULL, breaks=NULL, \dots, correction=NULL)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the multitype \eqn{J} function
    \eqn{J_{i\bullet}(r)}{Ji.(r)} will be computed.
    It must be a multitype point pattern (a marked point pattern
    whose marks are a factor). See under Details.
  }
  \item{i}{The type (mark value)
    of the points in \code{X} from which distances are measured.
    A character string (or something that will be converted to a
    character string).
    Defaults to the first level of \code{marks(X)}.
  }
  \item{eps}{A positive number.
    The resolution of the discrete approximation to Euclidean
    distance (see below). There is a sensible default.
  }
  \item{r}{numeric vector. The values of the argument \eqn{r}
    at which the function
    \eqn{J_{i\bullet}(r)}{Ji.(r)} should be evaluated.
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
    \code{"Hanisch"} and \code{"best"}.
  }
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).

  Essentially a data frame containing six numeric columns 
  \item{J}{the recommended
    estimator of \eqn{J_{i\bullet}(r)}{Ji.(r)},
    currently the Kaplan-Meier estimator.
  }
  \item{r}{the values of the argument \eqn{r} 
    at which the function \eqn{J_{i\bullet}(r)}{Ji.(r)} has been  estimated
  }
  \item{km}{the Kaplan-Meier 
    estimator of \eqn{J_{i\bullet}(r)}{Ji.(r)}
  }
  \item{rs}{the ``reduced sample'' or ``border correction''
    estimator of \eqn{J_{i\bullet}(r)}{Ji.(r)}
  }
  \item{han}{the Hanisch-style
    estimator of \eqn{J_{i\bullet}(r)}{Ji.(r)}
  }
  \item{un}{the ``uncorrected'' 
    estimator of \eqn{J_{i\bullet}(r)}{Ji.(r)}
    formed by taking the ratio of uncorrected empirical estimators
    of \eqn{1 - G_{i\bullet}(r)}{1 - Gi.(r)}
    and \eqn{1 - F_{\bullet}(r)}{1 - F.(r)}, see
    \code{\link{Gdot}} and \code{\link{Fest}}.
  }
  \item{theo}{the theoretical value of  \eqn{J_{i\bullet}(r)}{Ji.(r)}
    for a marked Poisson process, namely 1.
  }
  The result also has two attributes \code{"G"} and \code{"F"}
  which are respectively the outputs of \code{\link{Gdot}}
  and \code{\link{Fest}} for the point pattern.
}
\details{
  This function \code{Jdot} and its companions
  \code{\link{Jcross}} and \code{\link{Jmulti}}
  are generalisations of the function \code{\link{Jest}}
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
  The argument \code{i} will be interpreted as a
  level of the factor \code{X$marks}. (Warning: this means that
  an integer value \code{i=3} will be interpreted as the number 3,
  \bold{not} the 3rd smallest level.)
  
  The ``type \eqn{i} to any type'' multitype \eqn{J} function 
  of a stationary multitype point process \eqn{X}
  was introduced by Van lieshout and Baddeley (1999). It is defined by
  \deqn{J_{i\bullet}(r) = \frac{1 - G_{i\bullet}(r)}{1 -
      F_{\bullet}(r)}}{Ji.(r) = (1 - Gi.(r))/(1-F.(r))}
  where \eqn{G_{i\bullet}(r)}{Gi.(r)} is the distribution function of
  the distance from a type \eqn{i} point to the nearest other point
  of the pattern, and \eqn{F_{\bullet}(r)}{F.(r)} is the distribution
  function of the distance from a fixed point in space to the nearest
  point of the pattern.

  An estimate of \eqn{J_{i\bullet}(r)}{Ji.(r)}
  is a useful summary statistic in exploratory data analysis
  of a multitype point pattern. If the pattern is 
  a marked Poisson point process, then
  \eqn{J_{i\bullet}(r) \equiv 1}{Ji.(r) = 1}.
  If the subprocess of type \eqn{i} points is independent
  of the subprocess of points of all types not equal to \eqn{i},
  then \eqn{J_{i\bullet}(r)}{Ji.(r)} equals
  \eqn{J_{ii}(r)}{Jii(r)}, the ordinary \eqn{J} function
  (see \code{\link{Jest}} and Van Lieshout and Baddeley (1996))
  of the points of type \eqn{i}. 
  Hence deviations from zero of the empirical estimate of
  \eqn{J_{i\bullet} - J_{ii}}{Ji.-Jii} 
  may suggest dependence between types.

  This algorithm estimates \eqn{J_{i\bullet}(r)}{Ji.(r)} 
  from the point pattern \code{X}. It assumes that \code{X} can be treated
  as a realisation of a stationary (spatially homogeneous) 
  random spatial point process in the plane, observed through
  a bounded window.
  The window (which is specified in \code{X} as \code{Window(X)})
  may have arbitrary shape.
  Biases due to edge effects are
  treated in the same manner as in \code{\link{Jest}},
  using the Kaplan-Meier and border corrections.
  The main work is done by \code{\link{Gmulti}} and \code{\link{Fest}}.

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{J_{i\bullet}(r)}{Ji.(r)} should be evaluated. 
  The values of \eqn{r} must be increasing nonnegative numbers
  and the maximum \eqn{r} value must exceed the radius of the
  largest disc contained in the window.
}
\references{
  Van Lieshout, M.N.M. and Baddeley, A.J. (1996)
  A nonparametric measure of spatial interaction in point patterns.
  \emph{Statistica Neerlandica} \bold{50}, 344--361.

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
}
\seealso{
 \code{\link{Jcross}},
 \code{\link{Jest}},
 \code{\link{Jmulti}}
}
\examples{
     # Lansing woods data: 6 types of trees
   woods <- lansing

    \testonly{
        woods <- woods[seq(1,npoints(woods), by=30), ]
    }
    Jh. <- Jdot(woods, "hickory")
    plot(Jh.)
    # diagnostic plot for independence between hickories and other trees
    Jhh <- Jest(split(woods)$hickory)
    plot(Jhh, add=TRUE, legendpos="bottom")

    \dontrun{
    # synthetic example with two marks "a" and "b"
    pp <- runifpoint(30) \%mark\% factor(sample(c("a","b"), 30, replace=TRUE))
    J <- Jdot(pp, "a")
    }
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{nonparametric}

