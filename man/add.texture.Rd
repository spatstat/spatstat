\name{add.texture}
\alias{add.texture}
\title{
  Fill Plot With Texture
}
\description{
  Draws a simple texture inside a region on the plot.
}
\usage{
add.texture(W, texture = 4, spacing = NULL, ...)
}
\arguments{
  \item{W}{
    Window (object of class \code{"owin"}) inside which the
    texture should be drawn.
  }
  \item{texture}{
    Integer from 1 to 8 identifying the type of texture. See Details.
  }
  \item{spacing}{
    Spacing between elements of the texture, in units of the
    current plot.
  }
  \item{\dots}{
    Further arguments controlling the plot colour, line width etc.
  }
}
\details{
  The chosen texture, confined to the window \code{W},
  will be added to the current plot.
  The available textures are:
  \describe{
    \item{texture=1:}{
      Small crosses arranged in a square grid.
    }
    \item{texture=2:}{
      Parallel vertical lines.
    }
    \item{texture=3:}{
      Parallel horizontal lines.
    }
    \item{texture=4:}{ 
      Parallel diagonal lines at 45 degrees from the horizontal.
   }
    \item{texture=5:}{
      Parallel diagonal lines at 135 degrees from the horizontal.
    }
    \item{texture=6:}{
      Grid of horizontal and vertical lines.
    }
    \item{texture=7:}{
      Grid of diagonal lines at 45 and 135 degrees from the horizontal.
    }
    \item{texture=8:}{
      Grid of hexagons.
    }
  }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{owin}}, 
  \code{\link{plot.owin}},
  \code{\link{textureplot}},
  \code{\link{texturemap}}.
}
\examples{
  W <- Window(chorley)
  plot(W, main="")
  add.texture(W, 7)
}
\keyword{spatial}
\keyword{hplot}
