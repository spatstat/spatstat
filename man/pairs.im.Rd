\name{pairs.im}
\alias{pairs.im}
\title{
  Scatterplot Matrix for Pixel Images
}
\description{
  Produces a scatterplot matrix of the pixel values
  in two or more pixel images.
}
\usage{
\method{pairs}{im}(..., plot=TRUE)
}
\arguments{
  \item{\dots}{
    Any number of arguments, each of which is either
    a pixel image (object of class \code{"im"})
    or a named argument to be passed to \code{\link{pairs.default}}.
  }
  \item{plot}{
    Logical. If \code{TRUE}, the scatterplot matrix is plotted.
  }
}
\details{
  This is a method for the generic function \code{\link{pairs}}
  for the class of pixel images.
  
  It produces a square array of plot panels, in which each panel shows
  a scatterplot of the pixel values of one image against the
  corresponding pixel values of another image.

  At least two of the arguments \code{\dots}
  should be pixel images (objects of class \code{"im"}).
  Their spatial domains must overlap, but need not have the same
  pixel dimensions.

  First the pixel image domains are intersected, and converted to a
  common pixel resolution. Then the corresponding pixel values of each image are
  extracted. Then \code{\link{pairs.default}} is called to 
  plot the scatterplot matrix.

  Any arguments in \code{\dots} which are not pixel images will be
  passed to \code{\link{pairs.default}} to control the plot.
}
\section{Image or Contour Plots}{
  Since the scatterplots may show very dense concentrations of points,
  it may be useful to set \code{panel=panel.image}
  or \code{panel=panel.contour} to draw a colour image or contour plot
  of the kernel-smoothed density of the scatterplot in each panel.
  The argument \code{panel} is passed
  to \code{\link{pairs.default}}. See the help for 
  \code{\link{panel.image}} and \code{\link{panel.contour}}.
}
\section{Low Level Control of Graphics}{
  To control the appearance of the individual scatterplot panels,
  see \code{\link{pairs.default}}, \code{\link{points}}
  or \code{\link{par}}.
  To control the plotting symbol for the points in the scatterplot,
  use the arguments \code{pch}, \code{col}, \code{bg} as described
  under \code{\link{points}} (because 
  the default panel plotter is the function \code{\link{points}}).
  To suppress the tick marks on the plot axes,
  type \code{par(xaxt="n", yaxt="n")} before calling \code{pairs}.
}
\value{
  Invisible. A \code{data.frame} containing the
  corresponding pixel values for each image.
  The return value also belongs to the class \code{plotpairsim} which has
  a plot method, so that it can be re-plotted.
}
\seealso{
  \code{\link{pairs}},
  \code{\link{pairs.default}},
  \code{\link{panel.contour}},
  \code{\link{panel.image}},
  \code{\link{plot.im}},
  \code{\link{im}},
  \code{\link{par}}
}
\examples{
  X <- density(rpoispp(30))
  Y <- density(rpoispp(40))
  Z <- density(rpoispp(30))
  pairs(X,Y,Z)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{hplot}

