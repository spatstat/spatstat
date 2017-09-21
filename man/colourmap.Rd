\name{colourmap}
\alias{colourmap}
\title{Colour Lookup Tables}
\description{
  Create a colour map (colour lookup table).
}
\usage{
colourmap(col, ..., range=NULL, breaks=NULL, inputs=NULL)
}
\arguments{
  \item{col}{Vector of values specifying colours}
  \item{\dots}{Ignored.}
  \item{range}{
    Interval to be mapped.
    A numeric vector of length 2, specifying the endpoints of the
    range of values to be mapped.
    Incompatible with \code{breaks} or \code{inputs}.
  }
  \item{inputs}{
    Values to which the colours are associated.
    A factor or vector of the same length as \code{col}.
    Incompatible with \code{breaks} or \code{range}.
  }
  \item{breaks}{
    Breakpoints for the colour map.
    A numeric vector of length equal to \code{length(col)+1}.
    Incompatible with \code{range} or \code{inputs}.
  }
}
\details{
  A colour map is a mechanism for associating colours with data.
  It can be regarded as a function, mapping data to colours.

  The command \code{colourmap} creates an object representing
  a colour map, which can then be used to control the plot commands
  in the \pkg{spatstat} package. It can also be used to compute the
  colour assigned to any data value. 

  The argument \code{col} specifies the colours to which
  data values will be mapped. It should be a vector
  whose entries can be interpreted as colours by the standard
  \R graphics system. The entries can be string names of colours
  like \code{"red"}, or integers that refer to
  colours in the standard palette, or strings containing
  six-letter hexadecimal codes like \code{"#F0A0FF"}. 

  Exactly one of the arguments \code{range}, \code{inputs} or \code{breaks}
  must be specified by name.

  If \code{inputs} is given, then it should be a vector or factor,
  of the same length as \code{col}. The entries of \code{inputs} can be
  any atomic type (e.g. numeric, logical, character, complex) or factor
  values. The resulting colour map associates the value \code{inputs[i]}
  with the colour \code{col[i]}.

  If \code{range} is given, then it determines the interval of the real
  number line that will be mapped. It should be a numeric vector of
  length 2. 

  If \code{breaks} is given, then it determines the precise intervals
  of the real number line
  which are mapped to each colour. It should be a numeric vector,
  of length at least 2, with entries that are in increasing order.
  Infinite values are allowed. Any number in the range
  between \code{breaks[i]} and \code{breaks[i+1]} will be mapped to the
  colour \code{col[i]}. 

  The result is an object of class \code{"colourmap"}. 
  There are \code{print} and \code{plot} methods for this class.
  Some plot commands in the \pkg{spatstat} package accept an object
  of this class as a specification of the colour map.

  The result is also a function \code{f} which can be used to compute
  the colour assigned to any data value. 
  That is, \code{f(x)} returns the character value of the colour assigned
  to \code{x}. This also works for vectors of data values.
}
\value{
  A function, which is also an object of class \code{"colourmap"}.
}
\seealso{
  The plot method \code{\link{plot.colourmap}}.
  
  See the \R help file on
  \code{\link[grDevices:colors]{colours}} for information about the colours
  that \R recognises, and how to manipulate them.
  
  To make a smooth transition between colours, see
  \code{\link{interp.colourmap}}.
  To alter individual colour values, see
  \code{\link{tweak.colourmap}}.
    
  See \code{\link[spatstat:colourtools]{colourtools}}
  for more tools to manipulate colour values.

  See \code{\link{lut}} for lookup tables.
}
\examples{
  # colour map for real numbers, using breakpoints
  cr <- colourmap(c("red", "blue", "green"), breaks=c(0,5,10,15))
  cr
  cr(3.2)
  cr(c(3,5,7))
  # a large colour map
  co <- colourmap(rainbow(100), range=c(-1,1))
  co(0.2)
  # colour map for discrete set of values
  ct <- colourmap(c("red", "green"), inputs=c(FALSE, TRUE))
  ct(TRUE)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{color}
