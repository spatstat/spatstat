\name{beachcolours}
\alias{beachcolours}
\alias{beachcolourmap}
\title{
  Create Colour Scheme for a Range of Numbers
}
\description{
  Given a range of numerical values, this command creates a
  colour scheme that would be appropriate if the numbers were
  altitudes (elevation above or below sea level).
}
\usage{
beachcolours(range, sealevel = 0, monochrome = FALSE,
             ncolours = if (monochrome) 16 else 64,
             nbeach = 1)
beachcolourmap(range, ...) 
}
\arguments{
  \item{range}{
    Range of numerical values to be mapped.
    A numeric vector of length 2.
  }
  \item{sealevel}{
    Value that should be treated as zero.
    A single number,
    lying between \code{range[1]} and \code{range[2]}.
  }
  \item{monochrome}{
    Logical. If \code{TRUE} then a greyscale colour map is
    constructed. 
  }
  \item{ncolours}{
    Number of distinct colours to use.
  }
  \item{nbeach}{
    Number of colours that will be yellow.
  }
  \item{\dots}{Arguments passed to \code{beachcolours}.}
}
\details{
  Given a range of numerical values, these commands create a
  colour scheme that would be appropriate if the numbers were
  altitudes (elevation above or below sea level).

  Numerical values close to zero are portrayed in green (representing
  the waterline). Negative values are blue (representing water)
  and positive values are yellow to red (representing land).
  At least, these are the colours of land and sea in Western Australia.
  This colour scheme was proposed by Baddeley et al (2005).

  The function \code{beachcolours} returns these colours
  as a character vector, while \code{beachcolourmap}
  returns a colourmap object.

  The argument \code{range} should be a numeric vector of
  length 2 giving a range of numerical values.

  The argument \code{sealevel} specifies the height value that will
  be treated as zero, and mapped to the colour green.
  A vector of \code{ncolours} colours will be created,
  of which \code{nbeach} colours will be green.

  The argument \code{monochrome} is included
  for convenience when preparing publications.
  If \code{monochrome=TRUE} the colour map will be 
  a simple grey scale containing \code{ncolours}
  shades from black to white.
}
\value{
  For \code{beachcolours}, 
  a character vector of length \code{ncolours} specifying colour values.
  For \code{beachcolourmap}, a colour map (object of class \code{"colourmap"}).
}
\references{
  Baddeley, A., Turner, R., Moller, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.
}
\seealso{
  \code{\link{colourmap}}, 
  \code{\link[spatstat:colourtools]{colourtools}}.
}
\examples{
  plot(beachcolourmap(c(-2,2)))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{color}
