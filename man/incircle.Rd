\name{incircle}
\alias{incircle}
\title{Find Largest Circle Inside Window}
\description{
  Find the largest circle contained in a given window.
}
\usage{
incircle(W)
}
\arguments{
  \item{W}{A window (object of class \code{"owin"}).}
}
\details{
  Given a window \code{W} of any type and shape, 
  this function determines the largest circle
  that is contained inside \code{W}.

  For non-rectangular windows, the incircle is computed approximately
  by finding the maximum of the distance map (see
  \code{\link{distmap}}) of the complement of the window.
}
\value{
  A list with entries \code{x,y,r} giving the location \code{(x,y)}
  and radius \code{r} of the incircle.
}
\seealso{
  \code{\link{centroid.owin}}
}
\examples{
  W <- square(1)
  Wc <- incircle(W)
  plot(W)
  plot(disc(Wc$r, c(Wc$x, Wc$y)), add=TRUE)

  data(letterR)
  plot(letterR)
  Rc <- incircle(letterR)
  plot(disc(Rc$r, c(Rc$x, Rc$y)), add=TRUE)

  W <- as.mask(letterR)
  plot(W)
  Rc <- incircle(W)
  plot(disc(Rc$r, c(Rc$x, Rc$y)), add=TRUE)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}

