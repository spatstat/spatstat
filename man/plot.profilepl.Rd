\name{plot.profilepl}
\alias{plot.profilepl}
\title{
  Plot Profile Likelihood
}
\description{
  Plot the profile (pseudo) likelihood 
  against the irregular parameters, for a model that was fitted by
  maximum profile (pseudo)likelihood.
}
\usage{
\method{plot}{profilepl}(x, \dots, add = FALSE, main = NULL, tag = TRUE,
                         coeff = NULL, xvariable = NULL,
                         col = 1, lty = 1, lwd = 1,
                         col.opt = "green", lty.opt = 3, lwd.opt = 1)
}
\arguments{
  \item{x}{
    A point process model fitted by maximum profile (pseudo)likelihood.
    Object of class \code{"profilepl"}, obtained from \code{\link{profilepl}}.
  }
  \item{\dots}{
    Additional plot arguments passed to \code{\link[graphics]{plot.default}}
    and \code{\link[graphics]{lines}}.
  }
  \item{add}{
    Logical. If \code{TRUE}, the plot is drawn over the existing plot.
  }
  \item{main}{
    Optional. Main title for the plot. A character string or character vector.
  }
  \item{tag}{
    Logical value. If \code{TRUE} (the default), 
    when the plot contains multiple curves corresponding to different
    values of a parameter, each curve will be labelled 
    with the values of the irregular parameter.
  }
  \item{coeff}{
    Optional. If this is given, it should be a character string matching
    the name of one of the fitted model coefficients. This coefficient
    will then be plotted on the vertical axis.
  }
  \item{xvariable}{
    Optional. The name of the irregular parameter that should be plotted
    along the horizontal axis. The default is the first irregular
    parameter. 
  }
  \item{col,lty,lwd}{
    Graphical parameters (colour, line type, line width) for
    the curves on the plot. 
  }
  \item{col.opt, lty.opt, lwd.opt}{
    Graphical parameters for indicating the optimal parameter value.
  }
}
\details{
  This is the \code{\link[graphics]{plot}} method for the
  class \code{"profilepl"} of fitted point process models obtained
  by maximising the profile likelihood or profile pseudolikelihood.

  The default behaviour is to plot the profile likelihood or profile
  pseudolikelihood on the vertical axis, against the value of the
  irregular parameter on the horizontal axis.
  
  If there are several irregular parameters, then one of them
  is plotted on the horizontal axis, and the plot consists of many different
  curves, corresponding to different values of the other parameters.
  The parameter to be plotted on the horizontal axis is specified by the
  argument \code{xvariable}; the default is to use the parameter that was
  listed first in the original call to \code{\link{profilepl}}.

  If \code{coeff} is given, it should be the name of one of the
  fitted model coefficients \code{names(coef(as.ppm(x)))}. 
  The fitted value of that coefficient is plotted on the vertical axis.
}
\value{
  Null.
}
\references{
  Baddeley, A. and Turner, R. (2000)
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42}, 283--322.

  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  Chapman and Hall/CRC Press.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{profilepl}}
}
\examples{
  rstep <- if(interactive()) 0.005 else 0.02

  # one irregular parameter
  rr <- data.frame(r=seq(0.05,0.15, by=rstep))
  ps <- profilepl(rr, Strauss, cells)
  plot(ps)                      # profile pseudolikelihood 
  plot(ps, coeff="Interaction") # fitted interaction coefficient log(gamma)

  # two irregular parameters
  rs <- expand.grid(r=seq(0.05,0.15, by=rstep),sat=1:3)
  pg <- profilepl(rs, Geyer, cells)
  plot(pg) # profile pseudolikelihood against r for each value of 'sat'
  plot(pg, coeff="Interaction")
  plot(pg, xvariable="sat", col=ifelse(r < 0.1, "red", "green"))
}
\keyword{spatial}
\keyword{hplot}
\keyword{models}
