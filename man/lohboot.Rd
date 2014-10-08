\name{lohboot}
\alias{lohboot}
\title{Bootstrap Confidence Bands for Summary Function}
\description{
  Computes a bootstrap confidence band for a summary function
  of a point process.
}
\usage{
  lohboot(X,
          fun=c("pcf", "Kest", "Lest", "pcfinhom", "Kinhom", "Linhom"),
          \dots, nsim=200, confidence=0.95, global=FALSE, type=7)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{fun}{
    Name of the summary function for which confidence intervals are
    desired: one of the strings \code{"pcf"}, \code{"Kest"}, \code{"Lest"},
    \code{"pcfinhom"}, \code{"Kinhom"} or \code{"Linhom"}.
    Alternatively, the function itself; it must be
    one of the functions listed here.
  }
  \item{\dots}{
    Arguments passed to the corresponding local version of the summary
    function (see Details).
  }
  \item{nsim}{
    Number of bootstrap simulations.
  }
  \item{confidence}{
    Confidence level, as a fraction between 0 and 1.
  }
  \item{global}{
    Logical. If \code{FALSE} (the default), pointwise confidence intervals
    are constructed. If \code{TRUE}, a global (simultaneous) confidence band is
    constructed.
  }
  \item{type}{
    Integer. Argument passed to \code{\link[stats]{quantile}}
    controlling the way the quantiles are calculated.
  }
} 
\value{
  A function value table
  (object of class \code{"fv"})
  containing columns giving the estimate of the summary function,
  the upper and lower limits of the bootstrap confidence interval,
  and the theoretical value of the summary function for a Poisson process.
}
\details{
  This algorithm computes 
  confidence bands for the true value of the summary function
  \code{fun} using the bootstrap method of Loh (2008).

  If \code{fun="pcf"}, for example, the algorithm computes a pointwise
  \code{(100 * confidence)}\% confidence interval for the true value of
  the pair correlation function for the point process,
  normally estimated by \code{\link{pcf}}.
  It starts by computing the array of
  \emph{local} pair correlation functions,
  \code{\link{localpcf}}, of the data pattern \code{X}.
  This array consists of the contributions to the estimate of the
  pair correlation function from each
  data point. Then these contributions are resampled \code{nsim} times
  with replacement; from each resampled dataset the total contribution
  is computed, yielding \code{nsim} random pair correlation functions.
  The pointwise \code{alpha/2} and \code{1 - alpha/2} quantiles of
  these functions are computed, where \code{alpha = 1 - confidence}.
  The average of the local functions is also computed as an estimate
  of the pair correlation function.

  To control the estimation algorithm, use the 
  arguments \code{\dots}, which are passed to the local version
  of the summary function, as shown below:
  
  \tabular{ll}{
    \bold{fun} \tab \bold{local version} \cr
    \code{\link{pcf}} \tab \code{\link{localpcf}} \cr
    \code{\link{Kest}} \tab \code{\link{localK}} \cr
    \code{\link{Lest}} \tab \code{\link{localK}} \cr
    \code{\link{pcfinhom}} \tab \code{\link{localpcfinhom}} \cr
    \code{\link{Kinhom}} \tab \code{\link{localKinhom}} \cr
    \code{\link{Linhom}} \tab \code{\link{localKinhom}}
  }
  For \code{fun="Lest"}, the calculations are first performed
  as if \code{fun="Kest"}, and then the square-root transformation is
  applied to obtain the \eqn{L}-function. 

  Note that the confidence bands computed by 
  \code{lohboot(fun="pcf")} may not contain the estimate of the
  pair correlation function computed by \code{\link{pcf}},
  because of differences between the algorithm parameters
  (such as the choice of edge correction)
  in \code{\link{localpcf}} and \code{\link{pcf}}.
  If you are using \code{lohboot}, the
  appropriate point estimate of the pair correlation itself is
  the pointwise mean of the local estimates, which is provided
  in the result of \code{lohboot} and is shown in the default plot.

  If the confidence bands seem unbelievably narrow,
  this may occur because the point pattern has a hard core
  (the true pair correlation function is zero for certain values of
  distance) or because of an optical illusion when the
  function is steeply sloping (remember the width of the confidence
  bands should be measured \emph{vertically}).
  
  An alternative to \code{lohboot} is \code{\link{varblock}}.
}
\references{
  Loh, J.M. (2008)
  A valid and fast spatial bootstrap for correlation functions.
  \emph{The Astrophysical Journal}, \bold{681}, 726--734.
}
\seealso{
  Summary functions
  \code{\link{Kest}},
  \code{\link{pcf}},
  \code{\link{Kinhom}},
  \code{\link{pcfinhom}},
  \code{\link{localK}},
  \code{\link{localpcf}},
  \code{\link{localKinhom}},
  \code{\link{localpcfinhom}}.

  See \code{\link{varblock}} for an alternative bootstrap technique.
}
\examples{
  p <- lohboot(simdat, stoyan=0.5)
  plot(p)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
