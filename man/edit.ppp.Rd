\name{edit.ppp}
\alias{edit.ppp}
\alias{edit.im}
\title{
  Invoke Text Editor on Spatial Data
}
\description{
  Invokes a text editor allowing the user to inspect and change
  entries in a spatial dataset. 
}
\usage{
\method{edit}{ppp}(name, \dots)

\method{edit}{im}(name, \dots)
}
\arguments{
  \item{name}{
    Spatial dataset. An object of class \code{"ppp"} or \code{"im"}.
  }
  \item{\dots}{
    Other arguments passed to \code{\link[utils]{edit.data.frame}}.
  }
}
\details{
  The function \code{\link[utils]{edit}} is generic. These functions
  are methods for objects of class \code{"ppp"} and \code{"im"}.
  
  The spatial dataset \code{name} is converted to a data frame or array,
  and the text editor is invoked. The user can change the values of spatial
  coordinates or marks of the points in a point pattern,
  or the pixel values in an image.

  If \code{name} is a pixel image, it is converted to a matrix 
  and displayed in the same spatial orientation as if
  the image had been plotted.

  Note that the original object \code{name} is not changed;
  the function returns the edited dataset.
}
\value{
  Object of the same kind as \code{name} containing the edited data.
}
\author{
Adrian Baddeley
\email{Adrian.Baddeley@uwa.edu.au}
\url{http://www.maths.uwa.edu.au/~adrian/}
Rolf Turner
\email{r.turner@auckland.ac.nz}
and Ege Rubak
\email{rubak@math.aau.dk}
}
\seealso{
  \code{\link[utils]{edit.data.frame}}
}
\examples{
  if(interactive()) Z <- edit(cells)
}
\keyword{spatial}
\keyword{manip}

