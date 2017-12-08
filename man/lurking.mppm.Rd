\name{lurking.mppm}
\alias{lurking.mppm}
\title{
  Lurking Variable Plot for Multiple Point Patterns
}
\description{
  Generate a lurking variable plot of spatial point process residuals
  against a covariate, for a model fitted to several point patterns.
}
\usage{
\method{lurking}{mppm}(object, covariate, type="eem",
                       \dots,
                       separate = FALSE, 
                       plot.it = TRUE,
                       covname, oldstyle = FALSE, nx = 512, main="")
}
\arguments{
  \item{object}{
    The fitted model. An object of class \code{"mppm"} representing a
    point process model fitted to several point patterns.
  }
  \item{covariate}{
    The covariate to be used on the horizontal axis.
    Either an \code{expression} which can be evaluated in the
    original data, or a list of pixel images, one image for each
    point pattern in the original data.
  }
  \item{type}{
    String indicating the type of residuals or weights to be computed.
    Choices include \code{"eem"},
    \code{"raw"}, \code{"inverse"} and \code{"pearson"}.
    See \code{\link{diagnose.ppm}} for all possible choices.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{lurking.ppm}},
    including arguments controlling the plot.
  }
  \item{separate}{
    Logical value indicating whether to compute a separate lurking
    variable plot for each of the original point patterns.
    If \code{FALSE} (the default), a single lurking-variable plot
    is produced by combining residuals from all patterns.
  }
  \item{plot.it}{
    Logical value indicating whether 
    plots should be shown. If \code{plot.it=FALSE}, only
    the computed coordinates for the plots are returned.
    See \emph{Value}.
  }
  \item{covname}{
    A string name for the covariate, to be used in axis labels of plots.
  }
  \item{oldstyle}{
    Logical flag indicating whether error bounds should be plotted
    using the approximation given in the original paper
    (\code{oldstyle=TRUE}),
    or using the correct asymptotic formula (\code{oldstyle=FALSE}).
  }
  \item{nx}{
    Integer. Number of covariate values to be used in the plot.
  }
  \item{main}{
    Character string giving a main title for the plot.
  }
}
\details{
  This function generates a \sQuote{lurking variable} plot for a
  point process model fitted to several point patterns.
  Residuals from the model represented by \code{object}
  are plotted against the covariate specified by \code{covariate}.
  This plot can be used to reveal departures from the fitted model.

  The function \code{lurking} is generic. This is the method
  for the class \code{mppm}.
  The argument \code{object} must be a fitted point process model
  object of class \code{"mppm"}) produced by the
  model-fitting algorithm \code{\link{mppm}}.
}
\value{
  If \code{separate=FALSE} (the default), 
  the return value is an object 
  belonging to the class \code{"lurk"}, for which there
  are methods for \code{plot} and \code{print}.
  See \code{\link{lurking}} for details of the format.
  
  If \code{separate=TRUE}, the result is a list of such objects,
  and also belongs to the class \code{anylist} so that it can be
  printed and plotted. 
}
\author{
  \adrian,
  with thanks to Nicholas Read.
}
\seealso{
  \code{\link{lurking.ppm}}
}
\examples{
  fit <- mppm(Points ~ Image + Group, demohyper)
  lurking(fit, expression(Image), type="P")
  lurking(fit, expression(Image), type="P", separate=TRUE)
}
\keyword{spatial}
\keyword{models}
\keyword{hplot}
