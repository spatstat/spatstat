\name{colourtools}
\alias{colourtools} %DoNotExport
\alias{paletteindex}
\alias{rgb2hex}
\alias{col2hex}
\alias{paletteindex}
\alias{samecolour}
\title{
  Convert and Compare Colours in Different Formats
}
\description{
  These functions convert between different formats for specifying
  a colour in \R, and determine whether colours are equivalent.
}
\usage{
col2hex(x)
rgb2hex(v)
paletteindex(x)
samecolour(x,y)
}

\arguments{
  \item{x,y}{
    Any valid specification for a colour or sequence of colours
    accepted by \code{\link{col2rgb}}.
  }
  \item{v}{
    A numeric vector of length 3, giving the RGB values (0 to 255) of
    a single colour, or a 3-column matrix giving the RGB values of
    several colours.
  }
}
\details{
  \code{col2hex} converts colours specified in any format
  into their hexadecimal character codes.

  \code{rgb2hex} converts RGB colour values into their hexadecimal
  character codes.

  \code{paletteindex} checks whether the colour or colours specified
  by \code{x} are available in the default palette returned by
  \code{\link{palette}()}. If so, it returns the index or indices of
  the colours in the palette. If not, it returns \code{NA}.

  \code{samecolour} decides whether two colours \code{x} and \code{y}
  are equivalent.
}
\section{Warning}{
  \code{paletteindex("green")} returns \code{NA} because
  the green colour in the default palette is called \code{"green3"}.
}

\value{
  For \code{col2hex} and \code{rgb2hex}, a character vector containing
  hexadecimal colour codes.

  For \code{paletteindex}, an integer vector, possibly containing
  \code{NA} values.

  For \code{samecolour}, a logical value or logical vector.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{col2rgb}},
  \code{\link{palette}}

  See also the class of colour map objects:
  \code{\link{colourmap}}, 
  \code{\link{interp.colourmap}}, 
  \code{\link{tweak.colourmap}}.
}
\examples{
  samecolour("grey", "gray")
  paletteindex("grey")
}
\keyword{color}
