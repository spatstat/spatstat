\name{is.connected}  
\Rdversion{1.1}
\alias{is.connected}
\alias{is.connected.default}
\alias{is.connected.linnet}
\title{
  Determine Whether an Object is Connected
}
\description{
  Determine whether an object is 
  topologically connected.
}
\usage{
is.connected(X, \dots)

\method{is.connected}{default}(X, \dots)

\method{is.connected}{linnet}(X, \dots)
}
\arguments{
  \item{X}{
    A spatial object such as a pixel image (object of class \code{"im"}),
    a window (object of class \code{"owin"}) or a linear network
    (object of class \code{"linnet"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{connected}} to determine the
    connected components.
  }
}
\details{
  The command \code{is.connected(X)} returns \code{TRUE} if the object
  \code{X} consists of a single, topologically-connected piece,
  and returns \code{FALSE} if \code{X} consists of several pieces
  which are not joined together.
  
  The function \code{is.connected} is generic.
  The default method \code{is.connected.default} 
  works for many classes of objects, including windows (class \code{"owin"})
  and images (class \code{"im"}).
  There is a method for linear networks, \code{is.connected.linnet},
  described here, and a method for point patterns
  described in \code{\link{is.connected.ppp}}.
}
\value{
  A logical value.
}
\seealso{
  \code{\link{connected}},
  \code{\link{is.connected.ppp}}.
}
\examples{
  d <- distmap(cells, dimyx=256)
  X <- levelset(d, 0.07)
  plot(X)
  is.connected(X)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{math}
