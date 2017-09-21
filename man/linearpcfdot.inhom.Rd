\name{linearpcfdot.inhom}
\alias{linearpcfdot.inhom}
\title{
  Inhomogeneous Multitype
  Pair Correlation Function (Dot-type) for Linear Point Pattern
}
\description{
  For a multitype point pattern on a linear network,
  estimate the inhomogeneous multitype pair correlation function
  from points of type \eqn{i} to points of any type.
}
\usage{
linearpcfdot.inhom(X, i, lambdaI, lambdadot, r=NULL, \dots,
                   correction="Ang", normalise=TRUE)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the \eqn{i}-to-any pair correlation function
    \eqn{g_{i\bullet}(r)}{g[i.](r)} will be computed.
    An object of class \code{"lpp"} which 
    must be a multitype point pattern (a marked point pattern
    whose marks are a factor).
  }
  \item{i}{Number or character string identifying the type (mark value)
    of the points in \code{X} from which distances are measured.
    Defaults to the first level of \code{marks(X)}.
  }
  \item{lambdaI}{
    Intensity values for the points of type \code{i}. Either a numeric vector,
    a \code{function}, a pixel image
    (object of class \code{"im"} or \code{"linim"}) or
    a fitted point process model (object of class \code{"ppm"}
    or \code{"lppm"}).
  }
  \item{lambdadot}{
    Intensity values for all points of \code{X}. Either a numeric vector,
    a \code{function}, a pixel image
    (object of class \code{"im"} or \code{"linim"}) or
    a fitted point process model (object of class \code{"ppm"}
    or \code{"lppm"}).
  }
  \item{r}{numeric vector. The values of the argument \eqn{r}
    at which the function
    \eqn{g_{i\bullet}(r)}{g[i.](r)} should be evaluated.
    There is a sensible default.
    First-time users are strongly advised not to specify this argument.
    See below for important conditions on \eqn{r}.
  }
  \item{correction}{
    Geometry correction.
    Either \code{"none"} or \code{"Ang"}. See Details.
  }
  \item{\dots}{
    Arguments passed to \code{\link[stats]{density.default}}
    to control the kernel smoothing.
  }
  \item{normalise}{
    Logical. If \code{TRUE} (the default), the denominator of the estimator is 
    data-dependent (equal to the sum of the reciprocal intensities at
    the points of type \code{i}), which reduces the sampling variability.
    If \code{FALSE}, the denominator is the length of the network.
  }
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).
}
\details{
  This is a counterpart of the function \code{\link{pcfdot.inhom}} 
  for a point pattern on a linear network (object of class \code{"lpp"}).

  The argument \code{i} will be interpreted as
  levels of the factor \code{marks(X)}. 
  If \code{i} is missing, it defaults to the first
  level of the marks factor.

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{g_{i\bullet}(r)}{g[i.](r)}
  should be evaluated. 
  The values of \eqn{r} must be increasing nonnegative numbers
  and the maximum \eqn{r} value must not exceed the radius of the
  largest disc contained in the window.

  If \code{lambdaI} or \code{lambdadot} is a fitted point process model,
  the default behaviour is to update the model by re-fitting it to
  the data, before computing the fitted intensity.
  This can be disabled by setting \code{update=FALSE}.
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
 \code{\link{linearpcfcross.inhom}},
 \code{\link{linearpcfcross}},
 \code{\link{pcfcross.inhom}}.
}
\examples{
   lam <- table(marks(chicago))/(summary(chicago)$totlength)
   lamI <- function(x,y,const=lam[["assault"]]){ rep(const, length(x)) }
   lam. <- function(x,y,const=sum(lam)){ rep(const, length(x)) }

   g <- linearpcfdot.inhom(chicago, "assault", lamI, lam.)

   \dontrun{
     fit <- lppm(chicago, ~marks + x)
     linearpcfdot.inhom(chicago, "assault", fit, fit)
   }
}
\author{\adrian
  
  
}
\keyword{spatial}
\keyword{nonparametric}

