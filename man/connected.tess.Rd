\name{connected.tess}  
\Rdversion{1.1}
\alias{connected.tess}
\title{
  Connected Components of Tiles of a Tessellation
}
\description{
  Given a tessellation, find the topologically-connected pieces
  of each tile, and make a new tessellation using these pieces.
}
\usage{
\method{connected}{tess}(X, \dots)
}
\arguments{
  \item{X}{
    A tessellation (object of class \code{"tess"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} to determine the
    pixel resolution.
  }
}
\details{
  The function \code{connected} is generic.
  This function \code{connected.tess} is the method for tessellations.

  Given the tessellation \code{X}, the algorithm considers each tile
  of the tessellation, and identifies its connected components
  (topologically-connected pieces) using \code{\link{connected.owin}}.
  Each of these pieces is treated as a distinct tile
  and a new tessellation is made from these pieces.

  The result is another tessellation obtained by subdividing each tile
  of \code{X} into one or more new tiles.
}
\value{
  Another tessellation (object of class \code{"tess"}).
}
\seealso{
  \code{\link{connected.owin}}
}
\examples{
  BB <- grow.rectangle(Frame(letterR), 0.2)
  H <- tess(tiles=list(IN=letterR, OUT=complement.owin(letterR, BB)))
  opa <- par(mfrow=c(1,2))
  plot(H, do.col=TRUE)
  plot(connected(H), do.col=TRUE, col=2:4)
  par(opa)
}
\author{
  \adrian.
}
\keyword{spatial}
\keyword{math}
