\name{linearKcross}
\alias{linearKcross}
\title{
  Multitype K Function (Cross-type) for Linear Point Pattern
}
\description{
  For a multitype point pattern on a linear network,
  estimate the multitype \eqn{K} function
  which counts the expected number of points of type \eqn{j}
  within a given distance of a point of type \eqn{i}.
}
\usage{
linearKcross(X, i, j, r=NULL, \dots, correction="Ang")
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the cross type \eqn{K} function
    \eqn{K_{ij}(r)}{Kij(r)} will be computed.
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
    at which the \eqn{K}-function
    \eqn{K_{ij}(r)}{Kij(r)} should be evaluated.
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
  This is a counterpart of the function \code{\link{Kcross}} 
  for a point pattern on a linear network (object of class \code{"lpp"}).

  The arguments \code{i} and \code{j} will be interpreted as
  levels of the factor \code{marks(X)}. 
  If \code{i} and \code{j} are missing, they default to the first
  and second level of the marks factor, respectively.

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{K_{ij}(r)}{Kij(r)} should be evaluated. 
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
  The arguments \code{i} and \code{j} are interpreted as
  levels of the factor \code{marks(X)}. Beware of the usual
  trap with factors: numerical values are not
  interpreted in the same way as character values. 
}
\seealso{
 \code{\link{linearKdot}},
 \code{\link{linearK}}.
}
\examples{
   data(chicago)
   K <- linearKcross(chicago, "assault", "robbery")
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
}
\keyword{spatial}
\keyword{nonparametric}

