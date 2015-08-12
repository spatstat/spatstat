\name{linearKdot}
\alias{linearKdot}
\title{
  Multitype K Function (Dot-type) for Linear Point Pattern
}
\description{
  For a multitype point pattern on a linear network,
  estimate the multitype \eqn{K} function
  which counts the expected number of points (of any type)
  within a given distance of a point of type \eqn{i}.
}
\usage{
linearKdot(X, i, r=NULL, \dots, correction="Ang")
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the dot type \eqn{K} function
    \eqn{K_{i\bullet}(r)}{K[i.](r)} will be computed.
    An object of class \code{"lpp"} which 
    must be a multitype point pattern (a marked point pattern
    whose marks are a factor).
  }
  \item{i}{Number or character string identifying the type (mark value)
    of the points in \code{X} from which distances are measured.
    Defaults to the first level of \code{marks(X)}.
  }
  \item{r}{numeric vector. The values of the argument \eqn{r}
    at which the \eqn{K}-function
    \eqn{K_{i\bullet}(r)}{K[i.](r)} should be evaluated.
    There is a sensible default.
    First-time users are strongly advised not to specify this argument.
    See below for important conditions on \eqn{r}.
  }
  \item{correction}{
    Geometry correction.
    Either \code{"none"} or \code{"Ang"}. See Details.
  }
  \item{\dots}{Ignored.}
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).
}
\details{
  This is a counterpart of the function \code{\link{Kdot}} 
  for a point pattern on a linear network (object of class \code{"lpp"}).

  The argument \code{i} will be interpreted as
  levels of the factor \code{marks(X)}. 
  If \code{i} is missing, it defaults to the first
  level of the marks factor.

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{K_{i\bullet}(r)}{Ki.(r)} should be evaluated. 
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
 \code{\link{Kdot}},
 \code{\link{linearKcross}},
 \code{\link{linearK}}.
}
\examples{
   data(chicago)
   K <- linearKdot(chicago, "assault")
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{spatial}
\keyword{nonparametric}

