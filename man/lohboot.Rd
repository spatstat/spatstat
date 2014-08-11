\name{lohboot}
\alias{lohboot}
\title{Bootstrap Confidence Bands for Summary Function}
\description{
  Computes a bootstrap confidence band for a summary function
  of a point process.
}
\usage{
  lohboot(X,
          fun=c("pcf", "Kest", "pcfinhom", "Kinhom"),
          \dots, nsim=200, confidence=0.95, global=FALSE, type=7)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{fun}{
    Name of the summary function to be computed: one of the strings
    \code{"pcf"}, \code{"Kest"}, \code{"pcfinhom"} or \code{"Kinhom"}.
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
  confidence bands for the true value of the summary statistic
  \code{fun} using the bootstrap method of Loh (2008).

  If \code{fun="pcf"}, for example, the algorithm computes a pointwise
  \code{(100 * confidence)}\% confidence interval for the true value of
  the pair correlation function \code{\link{pcf}} for the point
  process. It starts by computing the array of
  \emph{local} pair correlation functions,
  \code{\link{localpcf}}, of the data pattern \code{X}.
  This array consists of the contributions to \code{\link{pcf}} from each
  data point. Then these contributions are resampled \code{nsim} times
  with replacement; from each resampled dataset the total contribution
  is computed, yielding \code{nsim} random pair correlation functions.
  The pointwise \code{alpha/2} and \code{1 - alpha/2} quantiles of
  these functions are computed, where \code{alpha = 1 - confidence}.

  To control the smoothing and estimation algorithm, use the 
  arguments \code{\dots}, which are passed to the local version
  of the summary function, as shown below:
  
  \tabular{ll}{
    \bold{fun} \tab \bold{local version} \cr
    \code{\link{pcf}} \tab \code{\link{localpcf}} \cr
    \code{\link{Kest}} \tab \code{\link{localK}} \cr
    \code{\link{pcfinhom}} \tab \code{\link{localpcfinhom}} \cr
    \code{\link{Kinhom}} \tab \code{\link{localKinhom}}
  }

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
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
