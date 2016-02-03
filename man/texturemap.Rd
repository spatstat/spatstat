\name{texturemap}
\alias{texturemap}
\title{
  Texture Map
}
\description{
  Create a map that associates data values with graphical textures.
}
\usage{
texturemap(inputs, textures, ...)
}
\arguments{
  \item{inputs}{
    A vector containing all the data values
    that will be mapped to textures.
  }
  \item{textures}{
    Optional. A vector of integer codes specifying the textures
    to which the \code{inputs} will be mapped.
  }
  \item{\dots}{
    Other graphics parameters such as \code{col}, \code{lwd}, \code{lty}.
  }
}
\details{
  A texture map is an association between data values and graphical
  textures. The command \code{texturemap} creates an object of class
  \code{"texturemap"} that represents a texture map.

  Once a texture map has been created, it can be applied to any
  suitable data to generate a texture plot of those data
  using \code{\link{textureplot}}.
  This makes it easy to ensure that
  the \emph{same} texture map is used in two different plots.
  The texture map can also be plotted in its own right.

  The argument \code{inputs} should be a vector containing all the
  possible data values (such as the levels of a factor) that are to be
  mapped.

  The \code{textures} should be integer values between 1 and 8,
  representing the eight possible textures
  described in the help for \code{\link{add.texture}}.
  The default is \code{textures = 1:n} where
  \code{n} is the length of \code{inputs}.
}
\value{
  An object of class \code{"texturemap"} representing the texture map.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}.
}
\seealso{
  \code{\link{textureplot}}
}
\examples{
   texturemap(letters[1:4], 2:5, col=1:4, lwd=2)
}
\keyword{spatial}
\keyword{hplot}
