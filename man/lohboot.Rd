\name{lohboot}
\alias{lohboot}
\title{Bootstrap Confidence Bands for Summary Function}
\description{
  Computes a bootstrap confidence band for a summary function
  of a point process.
}
\usage{
  lohboot(X,
          fun=c("pcf", "Kest", "Lest", "pcfinhom", "Kinhom", "Linhom",
                "Kcross", "Lcross", "Kdot", "Ldot",
                "Kcross.inhom", "Lcross.inhom"),
          \dots,
           block=FALSE, global=FALSE, basicboot=FALSE, Vcorrection=FALSE,
           confidence=0.95, nx = 4, ny = nx, nsim=200, type=7) 
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{fun}{
    Name of the summary function for which confidence intervals are
    desired: one of the strings \code{"pcf"}, \code{"Kest"}, \code{"Lest"},
    \code{"pcfinhom"}, \code{"Kinhom"} \code{"Linhom"}, \code{"Kcross"},
    \code{"Lcross"}, \code{"Kdot"}, \code{"Ldot"},
    \code{"Kcross.inhom"} or \code{"Lcross.inhom"}.
    Alternatively, the function itself; it must be
    one of the functions listed here.
  }
  \item{\dots}{
    Arguments passed to the corresponding local version of the summary
    function (see Details).
  }
  \item{block}{
    Logical value indicating whether to use Loh's block bootstrap
    as originally proposed. Default is \code{FALSE} for consistency
    with older code. See Details.
  }
  \item{global}{
    Logical. If \code{FALSE} (the default), pointwise confidence intervals
    are constructed. If \code{TRUE}, a global (simultaneous) confidence band is
    constructed.
  }
  \item{basicboot}{
    Logical value indicating whether to use the so-called basic
    bootstrap confidence interval. See Details.
  }
  \item{Vcorrection}{
    Logical value indicating whether to use a variance correction
    when \code{fun="Kest"} or \code{fun="Kinhom"}.
    See Details.
  }
  \item{confidence}{
    Confidence level, as a fraction between 0 and 1.
  }
  \item{nx,ny}{
    Integers.
    If \code{block=TRUE}, divide the window into \code{nx*ny} rectangles.
  }
  \item{nsim}{
    Number of bootstrap simulations.
  }
  \item{type}{
    Integer. Type of quantiles.
    Argument passed to \code{\link[stats]{quantile.default}}
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
  \code{fun} using the bootstrap method of Loh (2008)
  and a modification described in Baddeley, Rubak, Turner (2015).

  If \code{fun="pcf"}, for example, the algorithm computes a pointwise
  \code{(100 * confidence)}\% confidence interval for the true value of
  the pair correlation function for the point process,
  normally estimated by \code{\link{pcf}}.
  It starts by computing the array of
  \emph{local} pair correlation functions,
  \code{\link{localpcf}}, of the data pattern \code{X}.
  This array consists of the contributions to the estimate of the
  pair correlation function from each
  data point. 
  
  If  \code{block=FALSE}, these contributions are resampled \code{nsim} times
  with replacement as described in Baddeley, Rubak, Turner (2015);
  from each resampled dataset the total contribution
  is computed, yielding \code{nsim} random pair correlation functions.
  
  If \code{block=TRUE}, the calculation is performed as originally
  proposed by Loh (2008, 2010). The (bounding box of the) window is divided
  into \eqn{nx * ny} rectangles (blocks).
  The average contribution of a block
  is obtained by averaging the contribution of each point included in the block.
  Then, the average contributions on each block are resampled \code{nsim} times
  with replacement as described in Loh (2008) and Loh (2010);
  from each resampled dataset the total contribution
  is computed, yielding \code{nsim} random pair correlation functions.
  Notice that for non-rectangular windows any blocks not fully contained in the
  window are discarded before doing the resampling, so the effective number of
  blocks may be substantially smaller than \eqn{nx * ny} in this case.
  
  The pointwise \code{alpha/2} and \code{1 - alpha/2} quantiles of
  these functions are computed, where \code{alpha = 1 - confidence}.
  The average of the local functions is also computed as an estimate
  of the pair correlation function.

  There are several ways to define a bootstrap confidence interval. 
  If \code{basicbootstrap=TRUE},
  the so-called basic confidence bootstrap interval
  is used as described in Loh (2008).

  It has been noticed in Loh (2010) that
  when the intensity of the point process is unknown, 
  the bootstrap error estimate is larger than it should be.
  When the \eqn{K} function is used, 
  an adjustment procedure has been proposed in Loh (2010)
  that is used if \code{Vcorrection=TRUE}.
  In this case, the basic confidence bootstrap interval is implicitly used. 
  
  To control the estimation algorithm, use the 
  arguments \code{\dots}, which are passed to the local version
  of the summary function, as shown below:
  
  \tabular{ll}{
    \bold{fun} \tab \bold{local version} \cr
    \code{\link{pcf}} \tab \code{\link{localpcf}} \cr
    \code{\link{Kest}} \tab \code{\link{localK}} \cr
    \code{\link{Lest}} \tab \code{\link{localL}} \cr
    \code{\link{pcfinhom}} \tab \code{\link{localpcfinhom}} \cr
    \code{\link{Kinhom}} \tab \code{\link{localKinhom}} \cr
    \code{\link{Linhom}} \tab \code{\link{localLinhom}} \cr
    \code{\link{Kcross}} \tab \code{\link{localKcross}} \cr
    \code{\link{Lcross}} \tab \code{\link{localLcross}} \cr
    \code{\link{Kdot}} \tab \code{\link{localKdot}} \cr
    \code{\link{Ldot}} \tab \code{\link{localLdot}} \cr
    \code{\link{Kcross.inhom}} \tab \code{\link{localKcross.inhom}} \cr
    \code{\link{Lcross.inhom}} \tab \code{\link{localLcross.inhom}}
  }
  For \code{fun="Lest"}, the calculations are first performed
  as if \code{fun="Kest"}, and then the square-root transformation is
  applied to obtain the \eqn{L}-function. Similarly for
  \code{fun="Linhom", "Lcross", "Ldot", "Lcross.inhom"}.

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
  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  Chapman and Hall/CRC Press.

  Loh, J.M. (2008)
  A valid and fast spatial bootstrap for correlation functions.
  \emph{The Astrophysical Journal}, \bold{681}, 726--734.
  
  Loh, J.M. (2010)
  Bootstrapping an inhomogeneous point process.
  \emph{Journal of Statistical Planning and Inference}, \bold{140}, 734--749.
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
  \code{\link{localpcfinhom}},
  \code{\link{localKcross}},
  \code{\link{localKdot}},
  \code{\link{localLcross}},
  \code{\link{localLdot}}.
  \code{\link{localKcross.inhom}},
  \code{\link{localLcross.inhom}}.

  See \code{\link{varblock}} for an alternative bootstrap technique.
}
\examples{
  p <- lohboot(simdat, stoyan=0.5)
  g <- lohboot(simdat, stoyan=0.5, block=TRUE)
  g
  plot(g)
}
\author{
  \spatstatAuthors
  and Christophe Biscio.
}
\keyword{spatial}
\keyword{nonparametric}
