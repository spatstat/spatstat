\name{lut}
\alias{lut}
\title{Lookup Tables}
\description{
  Create a lookup table.
}
\usage{
lut(outputs, ..., range=NULL, breaks=NULL, inputs=NULL)
}
\arguments{
  \item{outputs}{Vector of output values}
  \item{\dots}{Ignored.}
  \item{range}{
    Interval of numbers to be mapped.
    A numeric vector of length 2, specifying the ends of the range of values
    to be mapped.
    Incompatible with \code{breaks} or \code{inputs}.
  }
  \item{inputs}{
    Input values to which the output values are associated.
    A factor or vector of the same length as \code{outputs}.
    Incompatible with \code{breaks} or \code{range}.
  }
  \item{breaks}{
    Breakpoints for the lookup table.
    A numeric vector of length equal to \code{length(outputs)+1}.
    Incompatible with \code{range} or \code{inputs}.
  }
}
\details{
  A lookup table is a function, mapping input values to output values.

  The command \code{lut} creates an object representing
  a lookup table, which can then be used to control various behaviour
  in the \pkg{spatstat} package. It can also be used to compute the
  output value assigned to any input value. 

  The argument \code{outputs} specifies the output values to which
  input data values will be mapped. It should be a vector of
  any atomic type (e.g. numeric, logical, character, complex) or factor
  values. 

  Exactly one of the arguments \code{range}, \code{inputs} or \code{breaks}
  must be specified by name.

  If \code{inputs} is given, then it should be a vector or factor,
  of the same length as \code{outputs}. The entries of \code{inputs} can be
  any atomic type (e.g. numeric, logical, character, complex) or factor
  values. The resulting lookup table associates the value \code{inputs[i]}
  with the value \code{outputs[i]}.

  If \code{range} is given, then it determines the interval of the real
  number line that will be mapped. It should be a numeric vector of
  length 2. 

  If \code{breaks} is given, then it determines intervals
  of the real number line
  which are mapped to each output value. It should be a numeric vector,
  of length at least 2, with entries that are in increasing order.
  Infinite values are allowed. Any number in the range
  between \code{breaks[i]} and \code{breaks[i+1]} will be mapped to the
  value \code{outputs[i]}. 

  The result is an object of class \code{"lut"}. 
  There is a \code{print} method for this class.
  Some plot commands in the \pkg{spatstat} package accept an object
  of this class as a specification of a lookup table.

  The result is also a function \code{f} which can be used to compute
  the output value assigned to any input data value. 
  That is, \code{f(x)} returns the output value assigned
  to \code{x}. This also works for vectors of input data values.
}
\value{
  A function, which is also an object of class \code{"lut"}.
}
\seealso{
  \code{\link{colourmap}}.
}
\examples{
  # lookup table for real numbers, using breakpoints
  cr <- lut(factor(c("low", "medium", "high")), breaks=c(0,5,10,15))
  cr
  cr(3.2)
  cr(c(3,5,7))
  # lookup table for discrete set of values
  ct <- lut(c(0,1), inputs=c(FALSE, TRUE))
  ct(TRUE)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{manip}
