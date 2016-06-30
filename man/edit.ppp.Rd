\name{edit.ppp}
\alias{edit.ppp}
\alias{edit.psp}
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

\method{edit}{psp}(name, \dots)

\method{edit}{im}(name, \dots)
}
\arguments{
  \item{name}{
    A spatial dataset
    (object of class \code{"ppp"}, \code{"psp"} or \code{"im"}).
  }
  \item{\dots}{
    Other arguments passed to \code{\link[utils]{edit.data.frame}}.
  }
}
\details{
  The function \code{\link[utils]{edit}} is generic. These functions
  are methods for spatial
  objects of class \code{"ppp"}, \code{"psp"} and \code{"im"}.

  The spatial dataset \code{name} is converted to a data frame or array,
  and the text editor is invoked. The user can change the values of spatial
  coordinates or marks of the points in a point pattern,
  or the coordinates or marks of the segments in a segment pattern,
  or the pixel values in an image. 
  The names of the columns of marks can also be edited.

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
\adrian


\rolf

and \ege

}
\seealso{
  \code{\link[utils]{edit.data.frame}},
  \code{\link{edit.hyperframe}}
}
\examples{
  if(interactive()) Z <- edit(cells)
}
\keyword{spatial}
\keyword{manip}

