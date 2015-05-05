\name{colourtools}
\alias{colourtools} %DoNotExport
\alias{paletteindex}
\alias{rgb2hex}
\alias{col2hex}
\alias{paletteindex}
\alias{samecolour}
\alias{complementarycolour}
\alias{interp.colours}
\alias{is.colour}
\alias{is.grey}
\alias{to.grey}
\title{
  Convert and Compare Colours in Different Formats
}
\description{
  These functions convert between different formats for specifying
  a colour in \R, determine whether colours are equivalent,
  and convert colour to greyscale.
}
\usage{
col2hex(x)
rgb2hex(v)
paletteindex(x)
samecolour(x,y)
complementarycolour(x)
interp.colours(x, length.out=512)
is.colour(x)
to.grey(x, weights=c(0.299, 0.587, 0.114))
is.grey(x)
}

\arguments{
  \item{x,y}{
    Any valid specification for a colour or sequence of colours
    accepted by \code{\link[grDevices]{col2rgb}}.
  }
  \item{v}{
    A numeric vector of length 3, giving the RGB values (0 to 255) of
    a single colour, or a 3-column matrix giving the RGB values of
    several colours.
  }
  \item{weights}{
    Numeric vector of length 3 giving 
    relative weights for the red, green, and blue
    channels respectively.
  }
  \item{length.out}{
    Integer. Length of desired sequence.
  }
}
\details{
  \code{is.colour(x)} can be applied to any kind of data \code{x}
  and returns \code{TRUE} if \code{x} can be interpreted as a colour or
  colours. The remaining functions expect data that can be interpreted
  as colours.
  
  \code{col2hex} converts colours specified in any format
  into their hexadecimal character codes.

  \code{rgb2hex} converts RGB colour values into their hexadecimal
  character codes.

  \code{paletteindex} checks whether the colour or colours specified
  by \code{x} are available in the default palette returned by
  \code{\link[grDevices]{palette}()}. If so, it returns the index or indices of
  the colours in the palette. If not, it returns \code{NA}.

  \code{samecolour} decides whether two colours \code{x} and \code{y}
  are equivalent.

  \code{is.grey} determines whether each entry of \code{x} is a
  greyscale colour, and returns a logical vector.
  
  \code{to.grey} converts the colour data in \code{x} to greyscale colours.
  Alternatively \code{x} can be an object of class \code{"colourmap"}
  and \code{to.grey(x)} is the modified colour map.

  \code{complementarycolour} replaces each colour by its
  complementary colour in RGB space (the colour obtained by replacing
  RGB values \code{(r, g, b)} by \code{(255-r, 255-g, 255-b)}).
  Alternatively \code{x} can be an object of class \code{"colourmap"}
  and \code{complementarycolour(x)} is the modified colour map.

  \code{interp.colours} interpolates between each successive pair of
  colours in a sequence of colours, to generate a more finely-spaced
  sequence. It uses linear interpolation in HSV space (with hue
  represented as a two-dimensional unit vector).
}
\section{Warning}{
  \code{paletteindex("green")} returns \code{NA} because
  the green colour in the default palette is called \code{"green3"}.
}

\value{
  For \code{col2hex}, \code{rgb2hex} and \code{to.grey},
  a character vector containing hexadecimal colour codes.

  For \code{paletteindex}, an integer vector, possibly containing
  \code{NA} values.

  For \code{samecolour} and \code{is.grey},
  a logical value or logical vector.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link[grDevices]{col2rgb}},
  \code{\link[grDevices]{palette}}

  See also the class of colour map objects:
  \code{\link{colourmap}}, 
  \code{\link{interp.colourmap}}, 
  \code{\link{tweak.colourmap}}.
}
\examples{
  samecolour("grey", "gray")
  paletteindex("grey")
  col2hex("orange")
  to.grey("orange")
  complementarycolour("orange")
  is.grey("lightgrey")
  is.grey(8)
  interp.colours(c("orange", "red", "violet"), 5)
}
\keyword{color}
