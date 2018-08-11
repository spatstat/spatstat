\name{symbolmap}
\alias{symbolmap}
\title{
  Graphics Symbol Map
}
\description{
  Create a graphics symbol map that associates data values with
  graphical symbols.
}
\usage{
symbolmap(\dots, range = NULL, inputs = NULL)
}
\arguments{
  \item{\dots}{
    Named arguments specifying the graphical parameters.
    See Details.
  }
  \item{range}{
    Optional. Range of numbers that are mapped.
    A numeric vector of length 2 giving the minimum and maximum
    values that will be mapped.
    Incompatible with \code{inputs}.
  }
  \item{inputs}{
    Optional. A vector containing all the data values
    that will be mapped to symbols.
    Incompatible with \code{range}.
  }
}
\details{
  A graphical symbol map is an association between
  data values and graphical symbols. 
  The command \code{symbolmap} creates an object of class
  \code{"symbolmap"} that represents a graphical symbol map.

  Once a symbol map has been created, it can be applied to any
  suitable data to generate a plot of those data.
  This makes it easy to ensure that
  the \emph{same} symbol map is used in two different plots.
  The symbol map can be plotted as a legend to the plots,
  and can also be plotted in its own right.
  
  The possible values of data that will be mapped
  are specified by \code{range} or \code{inputs}.
  \itemize{
    \item if \code{range} is given, it should be a numeric vector
    of length 2 giving the minimum and maximum values of the range
    of numbers that will be mapped. These limits must be finite.
    \item if \code{inputs} is given, it should be a vector
    of any atomic type (e.g. numeric, character, logical, factor).
    This vector contains all the possible data values
    that will be mapped.
    \item If neither \code{range} nor \code{inputs} is given,
    it is assumed that the possible values are real numbers.
  }
  The association of data values with graphical symbols
  is specified by the other arguments \code{\dots}
  which are given in \code{name=value} form.
  These arguments specify the kinds of symbols that will be
  used, the sizes of the symbols, and graphics parameters for
  drawing the symbols.

  Each graphics parameter can be either a single
  value, for example \code{shape="circles"},
  or a \code{function(x)} which determines the value
  of the graphics parameter as a function of the data \code{x},
  for example \code{shape=function(x) ifelse(x > 0, "circles", "squares")}.
  Colourmaps (see \code{\link{colourmap}}) are also acceptable
  because they are functions.
  
  Currently recognised graphics parameters, and their
  allowed values, are:
  \describe{
    \item{shape}{
      The shape of the symbol: currently
      either \code{"circles"}, \code{"squares"}, \code{"arrows"}
      or \code{NA}.
      This parameter takes precedence over \code{pch}.
    }
    \item{size}{
      The size of the symbol: a positive number or zero.
    }
    \item{pch}{
      Graphics character code:
      a positive integer, or a single character.
      See \code{\link[graphics]{par}}.
    }
    \item{cex}{
      Graphics character expansion factor.
    }
    \item{cols}{
      Colour of plotting characters.
    }
    \item{fg,bg}{
      Colour of foreground (or symbol border) and background
      (or symbol interior).
    }
    \item{col,lwd,lty}{
      Colour, width and style of lines.
    }
    \item{etch}{
      Logical. If \code{TRUE}, each symbol is surrounded
      by a border drawn in the opposite colour,
      which improves its visibility against the background.
      Default is \code{FALSE}.
    }
    \item{direction,headlength,headangle,arrowtype}{
      Numeric parameters of arrow symbols, applicable when
      \code{shape="arrows"}. Here \code{direction} is the direction
      of the arrow in degrees anticlockwise from the \eqn{x} axis;
      \code{headlength} is the length of the head of the arrow in
      coordinate units; \code{headangle} is the angle subtended by the point
      of the arrow; and \code{arrowtype} is an integer code
      specifying which ends of the shaft have arrowheads
      attached (0 means no arrowheads, 1 is an arrowhead at the start
      of the shaft, 2 is an arrowhead at the end of the shaft, and
      3 is arrowheads at both ends). 
    }
  }
  A vector of colour values is also acceptable for the arguments
  \code{col,cols,fg,bg} if
  \code{range} is specified.
}
\value{
  An object of class \code{"symbolmap"}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{plot.symbolmap}} to plot the symbol map itself.
  
  \code{\link{invoke.symbolmap}} to apply the symbol map to some data
  and plot the resulting symbols.
  
  \code{\link{update.symbolmap}} to change the symbol map.
}
\examples{
  g <- symbolmap(inputs=letters[1:10], pch=11:20)

  g1 <- symbolmap(range=c(0,100), size=function(x) x/50)

  g2 <- symbolmap(shape=function(x) ifelse(x > 0, "circles", "squares"),
                    size=function(x) sqrt(ifelse(x > 0, x/pi, -x)),
                    bg = function(x) ifelse(abs(x) < 1, "red", "black"))

  colmap <- colourmap(topo.colors(20), range=c(0,10))
  g3 <- symbolmap(pch=21, bg=colmap, range=c(0,10))
  plot(g3)
}
\keyword{spatial}
\keyword{hplot}
