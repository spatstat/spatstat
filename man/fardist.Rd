\name{fardist}
\alias{fardist}
\alias{fardist.ppp}
\alias{fardist.owin}
\title{
  Farthest Distance to Boundary of Window
}
\description{
  Computes the farthest distance from each pixel, or each data point,
  to the boundary of the window.
}
\usage{
  fardist(X, \dots)

  \method{fardist}{owin}(X, \dots, squared=FALSE)

  \method{fardist}{ppp}(X, \dots, squared=FALSE)
}
\arguments{
  \item{X}{
    A spatial object such as a window or point pattern.
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} to determine
    the pixel resolution, if required.
  }
  \item{squared}{
    Logical. If \code{TRUE}, the squared distances will be
    returned.
  }
}
\details{
  The function \code{fardist} is generic, with methods for
  the classes \code{owin} and \code{ppp}.

  For a window \code{W}, the command \code{fardist(W)}
  returns a pixel image in which the value at each pixel
  is the \emph{largest} distance from that pixel to the boundary of
  \code{W}.

  For a point pattern \code{X}, with window \code{W}, the
  command \code{fardist(X)} returns a numeric vector
  with one entry for each point of \code{X}, giving the
  largest distance from that data point to the boundary of \code{W}.
}
\value{
  For \code{fardist.owin}, a pixel image (object of class \code{"im"}).

  For \code{fardist.ppp}, a numeric vector.
}
\examples{
  fardist(cells)

  plot(FR <- fardist(letterR))
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\keyword{spatial}
\keyword{math}
