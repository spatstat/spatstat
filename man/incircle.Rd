\name{incircle}
\alias{incircle}
\alias{inradius}
\title{Find Largest Circle Inside Window}
\description{
  Find the largest circle contained in a given window.
}
\usage{
incircle(W)

inradius(W)
}
\arguments{
  \item{W}{A window (object of class \code{"owin"}).}
}
\details{
  Given a window \code{W} of any type and shape, 
  the function \code{incircle} determines the largest circle
  that is contained inside \code{W}, while \code{inradius} computes its
  radius only.

  For non-rectangular windows, the incircle is computed approximately
  by finding the maximum of the distance map (see
  \code{\link{distmap}}) of the complement of the window.
}
\value{
  The result of \code{incircle} is 
  a list with entries \code{x,y,r} giving the location \code{(x,y)}
  and radius \code{r} of the incircle.

  The result of \code{inradius} is the numerical value of radius.
}
\seealso{
  \code{\link{centroid.owin}}
}
\examples{
  W <- square(1)
  Wc <- incircle(W)
  plot(W)
  plot(disc(Wc$r, c(Wc$x, Wc$y)), add=TRUE)

  plot(letterR)
  Rc <- incircle(letterR)
  plot(disc(Rc$r, c(Rc$x, Rc$y)), add=TRUE)

  W <- as.mask(letterR)
  plot(W)
  Rc <- incircle(W)
  plot(disc(Rc$r, c(Rc$x, Rc$y)), add=TRUE)
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{math}

