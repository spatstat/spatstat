\name{Jcross}
\alias{Jcross}
\title{
  Multitype J Function (i-to-j)
}
\description{
  For a multitype point pattern, 
  estimate the multitype \eqn{J} function 
  summarising the interpoint dependence between
  points of type \eqn{i} and of type \eqn{j}.
}
\usage{
Jcross(X, i, j, eps=NULL, r=NULL, breaks=NULL, \dots, correction=NULL)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the multitype \eqn{J} function
    \eqn{J_{ij}(r)}{Jij(r)} will be computed.
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
  \item{eps}{A positive number.
    The resolution of the discrete approximation to Euclidean
    distance (see below). There is a sensible default.
  }
  \item{r}{Optional. Numeric vector. The values of the argument \eqn{r}
    at which the function
    \eqn{J_{ij}(r)}{Jij(r)} should be evaluated.
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
    Alternatively \code{correction="all"} selects all options.
  }
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).

  Essentially a data frame containing six numeric columns 
  \item{J}{the recommended
    estimator of \eqn{J_{ij}(r)}{Jij(r)},
    currently the Kaplan-Meier estimator.
  }
  \item{r}{the values of the argument \eqn{r} 
    at which the function \eqn{J_{ij}(r)}{Jij(r)} has been  estimated
  }
  \item{km}{the Kaplan-Meier 
    estimator of \eqn{J_{ij}(r)}{Jij(r)}
  }
  \item{rs}{the ``reduced sample'' or ``border correction''
    estimator of \eqn{J_{ij}(r)}{Jij(r)}
  }
  \item{han}{the Hanisch-style
    estimator of \eqn{J_{ij}(r)}{Jij(r)}
  }
  \item{un}{the ``uncorrected'' 
    estimator of \eqn{J_{ij}(r)}{Jij(r)}
    formed by taking the ratio of uncorrected empirical estimators
    of \eqn{1 - G_{ij}(r)}{1 - Gij(r)}
    and \eqn{1 - F_{j}(r)}{1 - Fj(r)}, see
    \code{\link{Gdot}} and \code{\link{Fest}}.
  }
  \item{theo}{the theoretical value of  \eqn{J_{ij}(r)}{Jij(r)}
    for a marked Poisson process, namely 1.
  }
  The result also has two attributes \code{"G"} and \code{"F"}
  which are respectively the outputs of \code{\link{Gcross}}
  and \code{\link{Fest}} for the point pattern.
}
\details{
  This function \code{Jcross} and its companions
  \code{\link{Jdot}} and \code{\link{Jmulti}}
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
  \bold{not} the 3rd smallest level).
  
  The ``type \eqn{i} to type \eqn{j}'' multitype \eqn{J} function 
  of a stationary multitype point process \eqn{X}
  was introduced by Van lieshout and Baddeley (1999). It is defined by
  \deqn{J_{ij}(r) = \frac{1 - G_{ij}(r)}{1 -
      F_{j}(r)}}{Jij(r) = (1 - Gij(r))/(1-Fj(r))}
  where \eqn{G_{ij}(r)}{Gij(r)} is the distribution function of
  the distance from a type \eqn{i} point to the nearest point of type \eqn{j},
  and \eqn{F_{j}(r)}{Fj(r)} is the distribution
  function of the distance from a fixed point in space to the nearest
  point of type \eqn{j} in the pattern.

  An estimate of \eqn{J_{ij}(r)}{Jij(r)}
  is a useful summary statistic in exploratory data analysis
  of a multitype point pattern. 
  If the subprocess of type \eqn{i} points is independent
  of the subprocess of points of type \eqn{j},
  then \eqn{J_{ij}(r) \equiv 1}{Jij(r) = 1}.
  Hence deviations of the empirical estimate of
  \eqn{J_{ij}}{Jij} from the value 1
  may suggest dependence between types.

  This algorithm estimates \eqn{J_{ij}(r)}{Jij(r)} 
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
  distance \eqn{r} at which \eqn{J_{ij}(r)}{Jij(r)} should be evaluated. 
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
  The arguments \code{i} and \code{j} are always interpreted as
  levels of the factor \code{X$marks}. They are converted to character
  strings if they are not already character strings.
  The value \code{i=1} does \bold{not}
  refer to the first level of the factor.
}
\seealso{
 \code{\link{Jdot}},
 \code{\link{Jest}},
 \code{\link{Jmulti}}
}
\examples{
     # Lansing woods data: 6 types of trees
    woods <- lansing
    \testonly{
       woods <- woods[seq(1,npoints(woods), by=30)]
    }
    Jhm <- Jcross(woods, "hickory", "maple")
    # diagnostic plot for independence between hickories and maples
    plot(Jhm)

    # synthetic example with two types "a" and "b"
    pp <- runifpoint(30) \%mark\% factor(sample(c("a","b"), 30, replace=TRUE))
    J <- Jcross(pp)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{nonparametric}

