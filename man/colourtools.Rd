\name{colourtools}
\alias{colourtools} %DoNotExport
\alias{paletteindex}
\alias{rgb2hex}
\alias{rgb2hsva}
\alias{col2hex}
\alias{paletteindex}
\alias{samecolour}
\alias{complementarycolour}
\alias{interp.colours}
\alias{is.colour}
\alias{is.grey}
\alias{to.grey}
\alias{to.opaque}
\alias{to.transparent}
\alias{to.saturated}
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
rgb2hex(v, maxColorValue=255)
rgb2hsva(red, green=NULL, blue=NULL, alpha=NULL, maxColorValue=255)
paletteindex(x)
samecolour(x,y)
complementarycolour(x)
interp.colours(x, length.out=512)
is.colour(x)
to.grey(x, weights=c(0.299, 0.587, 0.114), transparent=FALSE)
is.grey(x)
to.opaque(x)
to.transparent(x, fraction)
to.saturated(x, s=1)
}

\arguments{
  \item{x,y}{
    Any valid specification for a colour or sequence of colours
    accepted by \code{\link[grDevices]{col2rgb}}.
  }
  \item{v}{
    A numeric vector of length 3, giving the RGB values of
    a single colour, or a 3-column matrix giving the RGB values of
    several colours. Alternatively a vector of length 4 or a matrix with
    4 columns, giving the RGB and alpha (transparency) values.
  }
  \item{red,green,blue,alpha}{
    Arguments acceptable to \code{\link[grDevices]{rgb}}
    determining the red, green, blue channels and optionally the
    alpha (transparency) channel.
    Note that \code{red} can also be a matrix with 3 \bold{rows}
    giving the RGB values, or a matrix with 4 rows
    giving RGB and alpha values.
  }
  \item{maxColorValue}{
    Number giving the maximum possible value for the entries in
    \code{v} or \code{red,green,blue,alpha}.
  }
  \item{weights}{
    Numeric vector of length 3 giving 
    relative weights for the red, green, and blue
    channels respectively.
  }
  \item{transparent}{
    Logical value indicating whether transparent colours should
    be converted to transparent grey values (\code{transparent=TRUE})
    or converted to opaque grey values (\code{transparent=FALSE}, the
    default).
  }
  \item{fraction}{
    Transparency fraction. Numerical value or vector of values
    between 0 and 1, giving the opaqueness of a colour.
    A fully opaque colour has \code{fraction=1}.
  }
  \item{length.out}{
    Integer. Length of desired sequence.
  }
  \item{s}{
    Saturation value (between 0 and 1).
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
  character codes. It is a very minor extension to \code{\link[grDevices]{rgb}}.
  Arguments to \code{rgb2hex} should be similar to
  arguments to \code{\link[grDevices]{rgb}}.

  \code{rgb2hsva} converts RGB colour values into HSV colour values
  including the alpha (transparency) channel.
  It is an extension of \code{\link[grDevices]{rgb2hsv}}.
  Arguments to \code{rgb2hsva} should be similar to arguments to
  \code{\link[grDevices]{rgb2hsv}}.
  
  \code{paletteindex} checks whether the colour or colours specified
  by \code{x} are available in the default palette returned by
  \code{\link[grDevices]{palette}()}. If so, it returns the index or indices of
  the colours in the palette. If not, it returns \code{NA}.

  \code{samecolour} decides whether two colours \code{x} and \code{y}
  are equivalent.

  \code{is.grey} determines whether each entry of \code{x} is a
  greyscale colour, and returns a logical vector.
  
  \code{to.grey} converts the colour data in \code{x} to greyscale
  colours. Alternatively \code{x} can be an object of class \code{"colourmap"}
  and \code{to.grey(x)} is the modified colour map.

  \code{to.opaque} converts the colours in \code{x} to opaque
  (non-transparent) colours, and \code{to.transparent} converts them
  to transparent colours with a specified transparency value.
  Note that \code{to.transparent(x,1)} is equivalent to \code{to.opaque(x)}.

  For \code{to.grey}, \code{to.opaque} and
  \code{to.transparent}, if all the data in \code{x} specifies colours from the
  standard palette, and if the result would be equivalent to \code{x},
  then the result is identical to \code{x}.

  \code{to.saturated} converts each colour in \code{x} to
  its fully-saturated equivalent. For example, pink is mapped to red.
  Shades of grey are converted to black; white is unchanged.
  
  \code{complementarycolour} replaces each colour by its
  complementary colour in RGB space (the colour obtained by replacing
  RGB values \code{(r, g, b)} by \code{(255-r, 255-g, 255-b)}).
  The transparency value is not changed.
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
  For \code{col2hex} and \code{rgb2hex} 
  a character vector containing hexadecimal colour codes.

  For \code{to.grey}, \code{to.opaque} and \code{to.transparent},
  either a character vector containing hexadecimal colour codes,
  or a value identical to the input \code{x}.
  
  For \code{rgb2hsva}, a matrix with 3 or 4 rows containing
  HSV colour values.
  
  For \code{paletteindex}, an integer vector, possibly containing
  \code{NA} values.

  For \code{samecolour} and \code{is.grey},
  a logical value or logical vector.
}
\author{\adrian
  
  
  and \rolf
  
}
\seealso{
  \code{\link[grDevices]{col2rgb}},
  \code{\link[grDevices]{rgb2hsv}},
  \code{\link[grDevices]{palette}}.

  See also the class of colour map objects in the \pkg{spatstat} package:
  \code{\link{colourmap}}, 
  \code{\link{interp.colourmap}}, 
  \code{\link{tweak.colourmap}}.
}
\examples{
  samecolour("grey", "gray")
  paletteindex("grey")
  col2hex("orange")
  to.grey("orange")
  to.saturated("orange")
  complementarycolour("orange")
  is.grey("lightgrey")
  is.grey(8)
  to.transparent("orange", 0.5)
  to.opaque("red")
  interp.colours(c("orange", "red", "violet"), 5)
}
\keyword{color}
