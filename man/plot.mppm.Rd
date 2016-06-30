\name{plot.mppm}
\alias{plot.mppm}
\title{plot a Fitted Multiple Point Process Model}
\description{
  Given a point process model fitted to multiple point patterns
  by \code{\link{mppm}},
  compute spatial trend or conditional intensity surface of the model,
  in a form suitable for plotting, and (optionally) plot this
  surface.
}
\usage{
  \method{plot}{mppm}(x, \dots,
                trend=TRUE, cif=FALSE, se=FALSE,
                how=c("image", "contour", "persp"))
}
\arguments{
  \item{x}{
    A point process model fitted to multiple point patterns,
    typically obtained from
    the model-fitting algorithm \code{\link{mppm}}.
    An object of class \code{"mppm"}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.ppm}} or
    \code{\link{plot.anylist}} controlling the plot.
  }
  \item{trend}{
    Logical value indicating whether to plot the fitted trend.
  }
  \item{cif}{
    Logical value indicating whether to plot the fitted conditional intensity.
  }
  \item{se}{
    Logical value indicating whether to plot the standard error of the
    fitted trend.
  }
  \item{how}{
    Single character string indicating the style of plot to be performed. 
  }
}
\value{
  \code{NULL}.
}
\details{
  This is the \code{plot} method for the class \code{"mppm"}
  of point process models fitted to multiple point patterns
  (see \code{\link{mppm}}).

  It invokes \code{\link{subfits}} to compute the fitted model for
  each individual point pattern dataset, then calls
  \code{\link{plot.ppm}} to plot these individual models. These
  individual plots are displayed using \code{\link{plot.anylist}},
  which generates either a series of separate plot frames or an
  array of plot panels on a single page. 
}
\seealso{
  \code{\link{plot.ppm}},
  \code{\link{mppm}},
  \code{\link{plot.listof}}
}
\references{
  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  London: Chapman and Hall/CRC Press. 
}
\examples{
  # Synthetic data from known model
  n <- 9
  H <- hyperframe(V=1:n,
                  U=runif(n, min=-1, max=1))
  H$Z <- setcov(square(1))
  H$U <- with(H, as.im(U, as.rectangle(Z)))
  H$Y <- with(H, rpoispp(eval.im(exp(2+3*Z))))

  fit <- mppm(Y ~Z + U + V, data=H)

  plot(fit)
}
\author{
  \adrian, Ida-Maria Sintorn and Leanne Bischoff.
  Implemented by 
  \adrian
  
  
  \rolf
  
  and \ege
  
}
\keyword{spatial}
\keyword{hplot}
\keyword{models}
