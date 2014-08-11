\name{plot.layered}
\alias{plot.layered}
\title{
Layered Plot
}
\description{
  Generates a layered plot.
  The plot method for objects of class \code{"layered"}.
}
\usage{
\method{plot}{layered}(x, ..., which = NULL, plotargs = NULL,
                         add=FALSE, show.all=!add, main=NULL,
                         do.plot=TRUE)
}
\arguments{
  \item{x}{
    An object of class \code{"layered"}
    created by the function \code{\link{layered}}.
  }
  \item{\dots}{
    Arguments to be passed to the \code{plot} method
    for \emph{every} layer.
  }
  \item{which}{
    Subset index specifying which layers should be plotted. 
  }
  \item{plotargs}{
    Arguments to be passed to the  \code{plot} methods
    for individual layers. A list of lists of arguments of the
    form \code{name=value}.
  }
  \item{add}{Logical value indicating whether to add the graphics
    to an existing plot.
  }
  \item{show.all}{
     Logical value indicating whether the \emph{first} layer should be
     displayed in full (including the main title, bounding window,
     coordinate axes, colour ribbon, and so on). 
  }
  \item{main}{Main title for the plot}
  \item{do.plot}{Logical value indicating whether to actually do the plotting.}
}
\details{
  Layering is a simple mechanism for controlling
  a high-level plot that is composed of
  several successive plots, for example, a background and a foreground
  plot. The layering mechanism makes it easier to plot,
  to switch on or off the plotting of each individual layer,
  to control the plotting arguments that are passed to each layer,
  and to zoom in on a subregion.

  The layers of data to be plotted should first be converted
  into a single object of class \code{"layered"} using the
  function \code{\link{layered}}. Then the layers can be plotted using
  the method \code{plot.layered}.

  To zoom in on a subregion,
  apply the subset operator \code{\link{[.layered}} to \code{x}
  before plotting.

  Graphics parameters for each layer are determined by (in order of precedence)
  \code{\dots}, \code{plotargs}, and \code{\link{layerplotargs}(x)}.

  The argument \code{show.all} is recognised by many plot methods
  in \pkg{spatstat}. It determines whether a plot is drawn
  with all its additional components such as the main title, bounding window,
  coordinate axes, colour ribbons and legends. The default is
  \code{TRUE} for new plots and \code{FALSE} for added plots.
  
  In \code{plot.layered}, the argument \code{show.all} applies only to the
  \bold{first} layer. The subsequent layers are plotted
  with \code{show.all=FALSE}. 

  To override this, that is, if you really want to draw all
  the components of \bold{all} layers of \code{x},
  insert the argument \code{show.all=TRUE} in each
  entry of \code{plotargs} or \code{\link{layerplotargs}(x)}.
}
\value{
  (Invisibly) a list containing the return values from
  the plot commands for each layer. This list has
  an attribute \code{"bbox"} giving a bounding box for the entire plot.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{layered}},
  \code{\link{layerplotargs}},
  \code{\link{[.layered}},
  \code{\link{plot}}.
}
\examples{
   data(cells)
   D <- distmap(cells)
   L <- layered(D, cells)
   plot(L)
   plot(L, which = 2)
   plot(L, plotargs=list(list(ribbon=FALSE), list(pch=3, cols="white")))
   # plot a subregion
   plot(L[, square(0.5)])
}
\keyword{spatial}
\keyword{hplot}
