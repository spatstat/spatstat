\name{rQuasi}
\alias{rQuasi}
\title{
  Generate Quasirandom Point Pattern in Given Window
}
\description{
  Generates a quasirandom pattern of points in any two-dimensional window.
}
\usage{
rQuasi(n, W, type = c("Halton", "Hammersley"), ...)
}
\arguments{
  \item{n}{
    Maximum number of points to be generated.
  }
  \item{W}{
    Window (object of class \code{"owin"}) in which to generate
    the points.
  }
  \item{type}{
    String identifying the quasirandom generator.
  }
  \item{\dots}{
    Arguments passed to the quasirandom generator.
  }
}
\details{
  This function generates a quasirandom point pattern,
  using the quasirandom sequence generator
  \code{\link{Halton}} or \code{\link{Hammersley}}
  as specified.

  If \code{W} is a rectangle, 
  exactly \code{n} points will be generated.

  If \code{W} is not a rectangle, \code{n} points will be generated in the
  containing rectangle \code{as.rectangle(W)},
  and only the points lying inside \code{W} will be retained.
}
\value{
  Point pattern (object of class \code{"ppp"}) inside the window \code{W}.
}
\seealso{
  \code{\link{Halton}}
}
\examples{
   plot(rQuasi(256, letterR))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak.
}
\keyword{spatial}
\keyword{datagen}
