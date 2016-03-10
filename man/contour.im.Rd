\name{contour.im}
\alias{contour.im}
\title{Contour plot of pixel image}
\description{
  Generates a contour plot of a pixel image.
}
\usage{
   \method{contour}{im}(x, \dots, main,
            axes=FALSE, add=FALSE, col=par("fg"), 
            clipwin=NULL, show.all=!add, do.plot=TRUE)
}
\arguments{
  \item{x}{
    Pixel image to be plotted.
    An object of class \code{"im"}.
  }
  \item{main}{
    Character string to be displayed as the main title.
  }
  \item{axes}{
    Logical. If \code{TRUE}, coordinate axes are plotted
    (with tick marks) around a region slightly larger than the image window.
    If \code{FALSE} (the default), no axes are plotted,
    and a box is drawn tightly around the image window.
    Ignored if \code{add=TRUE}.
  }
  \item{add}{
    Logical. If \code{FALSE}, a new plot is created. If \code{TRUE},
    the contours are drawn over the existing plot.
  }
  \item{col}{
    Colour in which to draw the contour lines.
    Either a single value that can be interpreted as a colour value,
    or a \code{colourmap} object.
  }
  \item{clipwin}{
    Optional. A window (object of class \code{"owin"}).
    Only this subset of the data will be displayed.
  }
  \item{\dots}{
    Other arguments passed to \code{\link{contour.default}}
    controlling the contour plot; see Details.
  }
  \item{show.all}{
    Logical value indicating whether to display all plot elements
    including the main title, bounding box, and (if
    \code{axis=TRUE}) coordinate axis markings.
    Default is \code{TRUE} for new plots and \code{FALSE} for added plots.
  }
  \item{do.plot}{
    Logical value indicating whether to actually perform the plot.
  }
}
\details{
  This is a method for the generic \code{contour} function,
  for objects of the class \code{"im"}.
  
  An object of class \code{"im"} represents a pixel image;
  see \code{\link{im.object}}.

  This function displays the values of the pixel image \code{x}
  as a contour plot on the current plot device,
  using equal scales on the \eqn{x} and \eqn{y} axes.

  The appearance of the plot can be modified
  using any of the arguments listed in the help for
  \code{\link{contour.default}}.
  Useful ones include:
  \describe{
    \item{nlevels}{
      Number of contour levels to plot.
    }
    \item{drawlabels}{
      Whether to label the contour lines with text.
    }
    \item{col,lty,lwd}{
      Colour, type, and width of contour lines.
    }
  }
  See \code{\link{contour.default}} for a full list of these arguments.

  The defaults for any of the abovementioned arguments
  can be reset using \code{\link{spatstat.options}("par.contour")}.

  If \code{col} is a colour map (object of class \code{"colourmap"},
  see \code{\link{colourmap}}) then the contours will be plotted in
  different colours as determined by the colour map.
  The contour at level \code{z} will be plotted
  in the colour \code{col(z)} associated with this level in the colour map.
}
\value{
  none.
}
\examples{
   # an image
   Z <- setcov(owin())
   contour(Z, axes=TRUE)
   contour(Z)

   co <- colourmap(rainbow(100), range=c(0,1))
   contour(Z, col=co, lwd=2)
}
\seealso{
  \code{\link{im.object}},
  \code{\link{plot.im}},
  \code{\link{persp.im}}
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\keyword{spatial}
\keyword{hplot}

