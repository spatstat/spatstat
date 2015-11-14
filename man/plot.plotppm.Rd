\name{plot.plotppm}
\alias{plot.plotppm}
\title{Plot a plotppm Object Created by plot.ppm}
\description{
  The function plot.ppm produces objects which specify plots
  of fitted point process models.  The function plot.plotppm
  carries out the actual plotting of these objects.
}
\usage{
  \method{plot}{plotppm}(x, data = NULL, trend = TRUE, cif = TRUE,
             se = TRUE, pause = interactive(),
             how = c("persp", "image", "contour"),
             \dots, pppargs)
}
\arguments{
  \item{x}{
    An object of class \code{plotppm} produced by
    \code{\link{plot.ppm}()}
  }.
  \item{data}{
    The point pattern (an object of class \code{ppp})
    to which the point process model was fitted (by \code{\link{ppm}}).
  }
  \item{trend}{
    Logical scalar; should the trend component of
    the fitted model be plotted?
  }
  \item{cif}{
    Logical scalar; should the complete conditional
    intensity of the fitted model be plotted?
  }
  \item{se}{
    Logical scalar; should the estimated standard error
    of the fitted intensity be plotted?
  }
  \item{pause}{
    Logical scalar indicating whether to pause with a prompt 
    after each plot. Set \code{pause=FALSE} if plotting to a file.
  }
  \item{how}{
    Character string or character vector indicating the style or styles of
    plots to be performed.
  }
  \item{\dots}{
    Extra arguments to the plotting functions
    \code{\link{persp}}, \code{\link{image}} and \code{\link{contour}}.
  }
  \item{pppargs}{
    List of extra arguments passed to \code{\link{plot.ppp}}
    when displaying the original point pattern data.
  }
}
\details{
  If argument \code{data} is supplied then the point pattern will
  be superimposed on the image and contour plots.

  Sometimes a fitted model does not have a trend component, or the
  trend component may constitute all of the conditional intensity (if
  the model is Poisson).  In such cases the object \code{x} will not
  contain a trend component, or will contain only a trend component.
  This will also be the case if one of the arguments \code{trend}
  and \code{cif} was set equal to \code{FALSE} in the call to
  \code{plot.ppm()} which produced \code{x}.  If this is so then
  only the item which is present will be plotted.  Explicitly setting
  \code{trend=TRUE}, or \code{cif=TRUE}, respectively, will then give
  an error.
}
\value{
   None.
}

\section{Warning}{
  Arguments which are passed to \code{persp}, \code{image}, and
  \code{contour} via the \dots argument get passed to any of the
  other functions listed in the \code{how} argument, and won't be
  recognized by them.  This leads to a lot of annoying but harmless
  warning messages.  Arguments to \code{persp} may be supplied via
  \code{\link{spatstat.options}()} which alleviates the warning
  messages in this instance.
}

\author{\adrian
  
  
  and \rolf
  
}
\seealso{ \code{\link{plot.ppm}()} }
\examples{
 \dontrun{
 m <- ppm(cells ~ 1, Strauss(0.05))
 mpic <- plot(m)
 # Perspective plot only, with altered parameters:
  plot(mpic,how="persp", theta=-30,phi=40,d=4)
 # All plots, with altered parameters for perspective plot:
 op <- spatstat.options(par.persp=list(theta=-30,phi=40,d=4))
 plot(mpic)
 # Revert
 spatstat.options(op)
 }
}
\keyword{spatial}
\keyword{hplot}
\keyword{models}
