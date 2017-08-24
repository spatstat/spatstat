\name{plot.msr}
\alias{plot.msr}
\title{Plot a Signed or Vector-Valued Measure}
\description{
  Plot a signed measure or vector-valued measure.
}
\usage{
  \method{plot}{msr}(x, \dots,
                     add = FALSE,
                     how = c("image", "contour", "imagecontour"),
                     main = NULL,
                     do.plot = TRUE,
                     multiplot = TRUE,
                     massthresh = 0,
                     equal.markscale = FALSE,
                     equal.ribbon = FALSE)
}
\arguments{
  \item{x}{
    The signed or vector measure to be plotted.
    An object of class \code{"msr"} (see \code{\link{msr}}).
  }
  \item{\dots}{
    Extra arguments passed to \code{\link{Smooth.ppp}}
    to control the interpolation of the
    continuous density component of \code{x},
    or passed to \code{\link{plot.im}} or \code{\link{plot.ppp}}
    to control the appearance of the plot.
  }
  \item{add}{
    Logical flag; if \code{TRUE}, the graphics are added to the existing
    plot. If \code{FALSE} (the default) a new plot is initialised.
  }
  \item{how}{
    String indicating how to display the continuous density component.
  }
  \item{main}{
    String. Main title for the plot.
  }
  \item{do.plot}{
    Logical value determining whether to actually perform the plotting.
  }
  \item{multiplot}{
    Logical value indicating whether it is permissible to display
    a plot with multiple panels (representing different components of
    a vector-valued measure, or different types of points in a multitype
    measure.) 
  }
  \item{massthresh}{
    Threshold for plotting atoms.
    A single numeric value or \code{NULL}.
    If \code{massthresh=0} (the default) then only atoms with
    nonzero mass will be plotted.
    If \code{massthresh > 0} then only atoms whose absolute mass
    exceeds \code{massthresh} will be plotted.
    If \code{massthresh=NULL}, then all atoms of the measure will be plotted.
  }
  \item{equal.markscale}{
    Logical value indicating whether different panels should use the
    same symbol map (to represent the masses of atoms of the measure).
  }
  \item{equal.ribbon}{
    Logical value indicating whether different panels should use the
    same colour map (to represent the density values in the
    diffuse component of the measure).
  }
}
\value{
  (Invisible) colour map (object of class \code{"colourmap"}) for the
  colour image.
}
\details{
  This is the \code{plot} method for the class \code{"msr"}.
  
  The continuous density component of \code{x} is interpolated
  from the existing data by \code{\link{Smooth.ppp}},
  and then displayed as a colour image by \code{\link{plot.im}}.

  The discrete atomic component of \code{x} is then superimposed on this
  image by plotting the atoms as circles (for positive mass)
  or squares (for negative mass) by \code{\link{plot.ppp}}.
  By default, atoms with zero mass are not plotted at all. 

  To smooth both the discrete and continuous components,
  use   \code{\link{Smooth.msr}}.

  Use the argument \code{clipwin} to restrict the plot to a subset
  of the full data.

  To remove atoms with tiny masses, use the argument \code{massthresh}.
}
\seealso{
  \code{\link{msr}},
  \code{\link{Smooth.ppp}},
  \code{\link{Smooth.msr}},
  \code{\link{plot.im}},
  \code{\link{plot.ppp}}
}
\examples{
   X <- rpoispp(function(x,y) { exp(3+3*x) })
   fit <- ppm(X, ~x+y)
   rp <- residuals(fit, type="pearson")
   rs <- residuals(fit, type="score")

   plot(rp)
   plot(rs)
   plot(rs, how="contour")
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{hplot}
