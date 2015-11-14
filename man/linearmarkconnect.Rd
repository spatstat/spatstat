\name{linearmarkconnect}
\alias{linearmarkconnect}
\title{
  Mark Connection Function for Multitype Point Pattern on Linear Network
}
\description{
  For a multitype point pattern on a linear network,
  estimate the mark connection function
  from points of type \eqn{i} to points of type \eqn{j}.
}
\usage{
linearmarkconnect(X, i, j, r=NULL, \dots)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the mark connection function
    \eqn{p_{ij}(r)}{p[ij](r)} will be computed.
    An object of class \code{"lpp"} which 
    must be a multitype point pattern (a marked point pattern
    whose marks are a factor).
  }
  \item{i}{Number or character string identifying the type (mark value)
    of the points in \code{X} from which distances are measured.
    Defaults to the first level of \code{marks(X)}.
  }
  \item{j}{Number or character string identifying the type (mark value)
    of the points in \code{X} to which distances are measured.
    Defaults to the second level of \code{marks(X)}.
  }
  \item{r}{numeric vector. The values of the argument \eqn{r}
    at which the function
    \eqn{p_{ij}(r)}{p[ij](r)} should be evaluated.
    There is a sensible default.
    First-time users are strongly advised not to specify this argument.
    See below for important conditions on \eqn{r}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{linearpcfcross}}
    and \code{\link{linearpcf}}.
  }
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).
}
\details{
  This is a counterpart of the function \code{\link{markconnect}} 
  for a point pattern on a linear network (object of class \code{"lpp"}).

  The argument \code{i} will be interpreted as
  levels of the factor \code{marks(X)}. 
  If \code{i} is missing, it defaults to the first
  level of the marks factor.

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{p_{ij}(r)}{p[ij](r)}
  should be evaluated. 
  The values of \eqn{r} must be increasing nonnegative numbers
  and the maximum \eqn{r} value must not exceed the radius of the
  largest disc contained in the window.
}
\references{
  Baddeley, A, Jammalamadaka, A. and Nair, G. (to appear)
  Multitype point process analysis of spines on the
  dendrite network of a neuron.
  \emph{Applied Statistics} (Journal of the Royal Statistical
   Society, Series C), In press.
}
\section{Warnings}{
  The argument \code{i} is interpreted as a
  level of the factor \code{marks(X)}. Beware of the usual
  trap with factors: numerical values are not
  interpreted in the same way as character values. 
}
\seealso{
 \code{\link{linearpcfcross}},
 \code{\link{linearpcf}},
 \code{\link{linearmarkequal}},
 \code{\link{markconnect}}.
}
\examples{
   pab <- linearmarkconnect(chicago, "assault", "burglary")
  \dontrun{
   plot(alltypes(chicago, linearmarkconnect))
  }
}
\author{\adrian
  
  
}
\keyword{spatial}
\keyword{nonparametric}

