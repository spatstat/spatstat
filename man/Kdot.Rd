\name{Kdot}
\alias{Kdot}
\title{
  Multitype K Function (i-to-any)
}
\description{
  For a multitype point pattern, 
  estimate the multitype \eqn{K} function
  which counts the expected number of other points of the process
  within a given distance of a point of type \eqn{i}.
}
\usage{
Kdot(X, i, r=NULL, breaks=NULL, correction, ..., ratio=FALSE, from)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the multitype \eqn{K} function
    \eqn{K_{i\bullet}(r)}{Ki.(r)} will be computed.
    It must be a multitype point pattern (a marked point pattern
    whose marks are a factor). See under Details.
  }
  \item{i}{The type (mark value)
    of the points in \code{X} from which distances are measured.
    A character string (or something that will be converted to a
    character string).
    Defaults to the first level of \code{marks(X)}.
  }
  \item{r}{numeric vector. The values of the argument \eqn{r}
    at which the distribution function
    \eqn{K_{i\bullet}(r)}{Ki.(r)} should be evaluated.
    There is a sensible default.
    First-time users are strongly advised not to specify this argument.
    See below for important conditions on \eqn{r}.
  }
  \item{breaks}{
	This argument is for internal use only.
  }
  \item{correction}{
    A character vector containing any selection of the
    options \code{"border"}, \code{"bord.modif"},
    \code{"isotropic"}, \code{"Ripley"}, \code{"translate"},
    \code{"translation"},
    \code{"none"} or \code{"best"}.
    It specifies the edge correction(s) to be applied.
  }
  \item{\dots}{Ignored.}
  \item{ratio}{
    Logical. 
    If \code{TRUE}, the numerator and denominator of
    each edge-corrected estimate will also be saved,
    for use in analysing replicated point patterns.
  }
  \item{from}{An alternative way to specify \code{i}.}
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).

  Essentially a data frame containing numeric columns 
  \item{r}{the values of the argument \eqn{r} 
    at which the function \eqn{K_{i\bullet}(r)}{Ki.(r)} has been  estimated
  }
  \item{theo}{the theoretical value of  \eqn{K_{i\bullet}(r)}{Ki.(r)}
    for a marked Poisson process, namely \eqn{\pi r^2}{pi * r^2}
  }
  together with a column or columns named 
  \code{"border"}, \code{"bord.modif"},
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{K_{i\bullet}(r)}{Ki.(r)}
  obtained by the edge corrections named.

  If \code{ratio=TRUE} then the return value also has two
  attributes called \code{"numerator"} and \code{"denominator"}
  which are \code{"fv"} objects
  containing the numerators and denominators of each
  estimate of \eqn{K(r)}. 
}
\details{
  This function \code{Kdot} and its companions
  \code{\link{Kcross}} and \code{\link{Kmulti}}
  are generalisations of the function \code{\link{Kest}}
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
  level of the factor \code{X$marks}. 
  If \code{i} is missing, it defaults to the first
  level of the marks factor, \code{i = levels(X$marks)[1]}.
  
  The ``type \eqn{i} to any type'' multitype \eqn{K} function 
  of a stationary multitype point process \eqn{X} is defined so that
  \eqn{\lambda K_{i\bullet}(r)}{lambda Ki.(r)}
  equals the expected number of
  additional random points within a distance \eqn{r} of a
  typical point of type \eqn{i} in the process \eqn{X}.
  Here \eqn{\lambda}{lambda}
  is the intensity of the process,
  i.e. the expected number of points of \eqn{X} per unit area.
  The function \eqn{K_{i\bullet}}{Ki.} is determined by the 
  second order moment properties of \eqn{X}.

  An estimate of \eqn{K_{i\bullet}(r)}{Ki.(r)}
  is a useful summary statistic in exploratory data analysis
  of a multitype point pattern.
  If the subprocess of type \eqn{i} points were independent
  of the subprocess of points of all types not equal to \eqn{i},
  then \eqn{K_{i\bullet}(r)}{Ki.(r)} would equal \eqn{\pi r^2}{pi * r^2}.
  Deviations between the empirical \eqn{K_{i\bullet}}{Ki.} curve
  and the theoretical curve \eqn{\pi r^2}{pi * r^2} 
  may suggest dependence between types.

  This algorithm estimates the distribution function \eqn{K_{i\bullet}(r)}{Ki.(r)} 
  from the point pattern \code{X}. It assumes that \code{X} can be treated
  as a realisation of a stationary (spatially homogeneous) 
  random spatial point process in the plane, observed through
  a bounded window.
  The window (which is specified in \code{X} as \code{Window(X)})
  may have arbitrary shape.
  Biases due to edge effects are
  treated in the same manner as in \code{\link{Kest}},
  using the border correction.

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{K_{i\bullet}(r)}{Ki.(r)} should be evaluated. 
  The values of \eqn{r} must be increasing nonnegative numbers
  and the maximum \eqn{r} value must exceed the radius of the
  largest disc contained in the window.

  The pair correlation function can also be applied to the
  result of \code{Kdot}; see \code{\link{pcf}}.
}
\references{
  Cressie, N.A.C. \emph{Statistics for spatial data}.
    John Wiley and Sons, 1991.

  Diggle, P.J. \emph{Statistical analysis of spatial point patterns}.
    Academic Press, 1983.

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

}
\section{Warnings}{
  The argument \code{i} is interpreted as
  a level of the factor \code{X$marks}. It is converted to a character
  string if it is not already a character string.
  The value \code{i=1} does \bold{not}
  refer to the first level of the factor.

  The reduced sample estimator of \eqn{K_{i\bullet}}{Ki.} is pointwise approximately 
  unbiased, but need not be a valid distribution function; it may 
  not be a nondecreasing function of \eqn{r}. Its range is always 
  within \eqn{[0,1]}.
}
\seealso{
 \code{\link{Kdot}},
 \code{\link{Kest}},
 \code{\link{Kmulti}},
 \code{\link{pcf}}
}
\examples{
     # Lansing woods data: 6 types of trees
    woods <- lansing

    \testonly{woods <- woods[seq(1, npoints(woods), by=80)]}

    Kh. <- Kdot(woods, "hickory") 
    # diagnostic plot for independence between hickories and other trees
    plot(Kh.)

    \dontrun{
    # synthetic example with two marks "a" and "b"
    pp <- runifpoispp(50)
    pp <- pp \%mark\% factor(sample(c("a","b"), npoints(pp), replace=TRUE))
    K <- Kdot(pp, "a")
    }
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{nonparametric}

