\name{Window}
\alias{Window}
\alias{Window<-}
\alias{Window.ppp}
\alias{Window<-.ppp}
\alias{Window.psp}
\alias{Window<-.psp}
\alias{Window.im}
\alias{Window<-.im}
\title{
  Extract or Change the Window of a Spatial Object
}
\description{
  Given a spatial object (such as a point pattern or pixel image)
  in two dimensions, these functions extract or change the window
  in which the object is defined.
}
\usage{
   Window(X, \dots)

   Window(X, \dots) <- value

   \method{Window}{ppp}(X, \dots)

   \method{Window}{ppp}(X, \dots) <- value

   \method{Window}{psp}(X, \dots)

   \method{Window}{psp}(X, \dots) <- value

   \method{Window}{im}(X, \dots)

   \method{Window}{im}(X, \dots) <- value
}
\arguments{
  \item{X}{
    A spatial object such as a point pattern, line segment pattern
    or pixel image.
  }
  \item{\dots}{
    Extra arguments. They are ignored by all the methods listed here.
  }
  \item{value}{
    Another window (object of class \code{"owin"}) to be used as the
    window for \code{X}.
  }
}
\details{
  The functions \code{Window} and \code{Window<-} are generic.

  \code{Window(X)} extracts the spatial window in which \code{X} is
  defined.

  \code{Window(X) <- W} changes the window in which \code{X} is defined
  to the new window \code{W}, and \emph{discards any data outside} \code{W}.
  In particular:
  \itemize{
    \item If \code{X} is a point pattern (object of class \code{"ppp"})
    then \code{Window(X) <- W} discards any points of \code{X} which
    fall outside \code{W}.
    \item If \code{X} is a line segment pattern (object of class
    \code{"psp"}) then \code{Window(X) <- W} clips the segments of \code{X}
    to the boundaries of \code{W}.
    \item If \code{X} is a pixel image (object of class \code{"im"})
    then \code{Window(X) <- W} has the effect that pixels
    lying outside \code{W} are retained but their pixel values
    are set to \code{NA}.
  }

  Many other classes of spatial object have a method
  for \code{Window}, but not \code{Window<-}.
  See \code{\link{Window.ppm}}.
}
\value{
  The result of \code{Window} is a window (object of class
  \code{"owin"}).

  The result of \code{Window<-} is the updated object \code{X},
  of the same class as \code{X}.
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\seealso{
  \code{\link{Window.ppm}}
}
\examples{
   ## point patterns
   Window(cells)
   X <- demopat
   Window(X)
   Window(X) <- as.rectangle(Window(X))

   ## line segment patterns
   X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
   Window(X)
   Window(X) <- square(0.5)

   ## images
   Z <- setcov(owin())
   Window(Z)
   Window(Z) <- square(0.5)
}
\keyword{spatial}
\keyword{manip}
