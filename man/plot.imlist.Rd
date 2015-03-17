\name{plot.imlist}
\alias{plot.imlist}
\alias{image.imlist}
\alias{image.listof}
\title{Plot a List of Images}
\description{
  Plots an array of pixel images.
}
\usage{
  \method{plot}{imlist}(x, \dots, plotcommand="image",
                          equal.ribbon=FALSE, ribmar=NULL)

  \method{image}{imlist}(x, \dots, equal.ribbon=FALSE, ribmar=NULL)

  \method{image}{listof}(x, \dots, equal.ribbon=FALSE, ribmar=NULL)
}
\arguments{
  \item{x}{
    An object of the class \code{"imlist"}
    representing a list of pixel images.
    Alternatively \code{x} may belong to the outdated class \code{"listof"}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.solist}} to control the
    spatial arrangement of panels, and arguments passed to \code{\link{plot.im}}
    to control the display of each panel.
  }
  \item{equal.ribbon}{
    Logical. If \code{TRUE}, the colour maps of all the images will be
    the same. If \code{FALSE}, the colour map of each
    image is adjusted to the range of values of that image.
  }
  \item{ribmar}{
    Numeric vector of length 4 specifying the margins around the
    colour ribbon, if \code{equal.ribbon=TRUE}. Entries in the vector
    give the margin at the bottom, left, top, and right respectively,
    as a multiple of the height of a line of text.
  }
  \item{plotcommand}{
    Character string giving the name of a function
    to be used to display each image.
    Recognised by \code{plot.imlist} only.
  }
}
\value{
  Null.
}
\details{
  These are methods for the generic plot commands
  \code{plot} and \code{image} for the class \code{"imlist"}.
  They are currently identical.
  
  An object of class \code{"imlist"} represents a list of pixel images.
  (The outdated class \code{"listof"} is also handled.)

  Each entry in the list \code{x} will be displayed as a pixel image,
  in an array of panels laid out on the same graphics display,
  using \code{\link{plot.solist}}. Individual panels are plotted by
  \code{\link{plot.im}}.

  If \code{equal.ribbon=FALSE} (the default), 
  the images are rendered using different colour maps,
  which are displayed as colour ribbons beside each image.
  If \code{equal.ribbon=TRUE}, the images are rendered using the
  same colour map, and a single colour ribbon will be displayed at the right
  side of the array. 
  The colour maps and the placement of the colour ribbons are
  controlled by arguments \code{\dots} passed to \code{\link{plot.im}}.
}
\seealso{
  \code{\link{plot.solist}},
  \code{\link{plot.im}}
}
\examples{
 D <- density(split(amacrine))
 image(D, equal.ribbon=TRUE, main="", col.ticks="red", col.axis="red")
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{hplot}
