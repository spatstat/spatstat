\name{plot.ppm}
\alias{plot.ppm}
\title{plot a Fitted Point Process Model}
\description{
  Given a fitted point process model obtained by \code{\link{ppm}},
  create spatial trend and conditional intensity surfaces of the model,
  in a form suitable for plotting, and (optionally) plot these
  surfaces.
}
\usage{
  \method{plot}{ppm}(x, ngrid = c(40,40), superimpose = TRUE, 
                  trend = TRUE, cif = TRUE, se = TRUE, pause = interactive(),
                  how=c("persp","image", "contour"), plot.it = TRUE,
                  locations = NULL, covariates=NULL, \dots)
}
\arguments{
  \item{x}{
    A fitted point process model, typically obtained from
    the model-fitting algorithm \code{\link{ppm}}.
    An object of class \code{"ppm"}.
  }
  \item{ngrid}{
    The dimensions for a grid on which to evaluate,
    for plotting, the spatial trend and conditional intensity.
    A vector of 1 or 2 integers. If it is of length 1,
    \code{ngrid} is replaced by \code{c(ngrid,ngrid)}.
  }
  \item{superimpose}{
    logical flag; if \code{TRUE} (and if \code{plot=TRUE}) the
    original data point pattern will be superimposed on the plots.
  }
  \item{trend}{
    logical flag; if \code{TRUE}, the spatial trend surface will be produced.
  }
  \item{cif}{
    logical flag; if \code{TRUE}, the conditional intensity surface will be
    produced.
  }
  \item{se}{
    logical flag; if \code{TRUE}, the estimated standard error of the
    spatial trend surface will be produced.
  }
  \item{pause}{
    logical flag indicating whether to pause with a prompt 
    after each plot. Set \code{pause=FALSE} if plotting to a file.
    (This flag is ignored if \code{plot=FALSE}).
  }
  \item{how}{
    character string or character vector indicating the style or styles of
    plots to be performed. Ignored if \code{plot=FALSE}.  
  }
  \item{plot.it}{
    logical scalar; should a plot be produced immediately?
  }
  \item{locations}{
    If present, this determines the locations of the pixels
    at which predictions are computed. It must be a binary pixel image
    (an object of class \code{"owin"} with type \code{"mask"}).
    (Incompatible with \code{ngrid}).
  }
  \item{covariates}{
    Values of external covariates required by the fitted model.
    Passed to \code{\link{predict.ppm}}.
  }
  \item{\dots}{
    extra arguments to the plotting functions \code{\link{persp}},
    \code{\link{image}} and \code{\link{contour}}.
  }
}
\value{
  An object of class \code{plotppm}.  Such objects may be plotted by
  \code{\link{plot.plotppm}()}.

  This is a list with components named \code{trend} and \code{cif},
  either of which may
  be missing. They will be missing if the corresponding component
  does not make sense for the model, or if the corresponding
  argument was set equal to \code{FALSE}.

  Both \code{trend} and \code{cif} are lists of images.
  If the model is an unmarked point process, then they are lists of
  length 1, so that \code{trend[[1]]} is an image of the spatial trend
  and \code{cif[[1]]} is an image of the conditional intensity.

  If the model is a marked point process, then \code{trend[[i]]}
  is an image of the spatial trend for the mark \code{m[i]},
  and \code{cif[[1]]} is an image of the conditional intensity
  for the mark \code{m[i]}, where \code{m} is the vector of levels
  of the marks.
}
\details{
  This is the \code{plot} method for the class \code{"ppm"}
  (see \code{\link{ppm.object}} for details of this class).
  
  It invokes \code{\link{predict.ppm}} to compute the spatial 
  trend and conditional intensity of the fitted point process model.
  See \code{\link{predict.ppm}} for more explanation about spatial trend
  and conditional intensity.
 
  The default action is to create a rectangular grid
  of points in (the bounding box of) the observation window of
  the data point pattern, and evaluate the spatial trend and
  conditional intensity of the fitted spatial point process model
  \code{x} at these locations.  If the argument \code{locations=}
  is supplied, then the spatial trend
  and conditional intensity are calculated at the grid of points
  specified by this argument.

  The argument \code{locations}, if present, should be a
  binary image mask (an object of class \code{"owin"}
  and type \code{"mask"}). This determines a rectangular grid
  of locations, or a subset of such a grid, at which predictions
  will be computed. Binary image masks
  are conveniently created using \code{\link{as.mask}}.

  The argument \code{covariates} gives the values of any spatial covariates
  at the prediction locations.
  If the trend formula in the fitted model 
  involves spatial covariates (other than
  the Cartesian coordinates \code{x}, \code{y})
  then \code{covariates} is required.

  The argument \code{covariates} has the same format and interpretation
  as in \code{\link{predict.ppm}}. It may be
  either a data frame (the number of whose rows must match
  the number of pixels in \code{locations} multiplied by the number of
  possible marks in the point pattern),  or a list of images.
  If argument \code{locations}
  is not supplied, and \code{covariates} \bold{is} supplied, then
  it \bold{must} be a list of images.

  If the fitted model was a marked (multitype) point process, then 
  predictions are made for each possible mark value in turn.
 
  If the fitted model had no spatial trend, then the default is 
  to omit calculating this (flat) surface, unless \code{trend=TRUE}
  is set explicitly.
 
  If the fitted model was Poisson, so that there were no spatial interactions,
  then the conditional intensity and spatial trend are identical, and the
  default is to omit the conditional intensity, unless \code{cif=TRUE} is set
  explicitly.

  If \code{plot.it=TRUE} then \code{\link{plot.plotppm}()} is called
  upon to plot the class \code{plotppm} object which is produced.
  (That object is also returned, silently.)
  
  Plots are produced successively using \code{\link{persp}},
  \code{\link{image}} and \code{\link{contour}} (or only a
  selection of these three, if \code{how} is given).  Extra
  graphical parameters controlling the display may be passed
  directly via the arguments \code{...} or indirectly reset using
  \code{\link{spatstat.options}}.
}
\seealso{
  \code{\link{plot.plotppm}},
  \code{\link{ppm}},
  \code{\link{ppm.object}},
  \code{\link{predict.ppm}},
  \code{\link{print.ppm}},
  \code{\link{persp}},
  \code{\link{image}},
  \code{\link{contour}},
  \code{\link{plot}},
  \code{\link{spatstat.options}}
}
\section{Warnings}{
  See warnings in \code{\link{predict.ppm}}.
}
\examples{
 \dontrun{
 data(cells)
 m <- ppm(cells, ~1, Strauss(0.05))
 pm <- plot(m) # The object ``pm'' will be plotted as well as saved
               # for future plotting.
 }
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{hplot}
\keyword{models}
