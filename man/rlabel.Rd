\name{rlabel}
\alias{rlabel}
\title{Random Re-Labelling of Point Pattern}
\description{
  Randomly allocates marks to a point pattern,
  or permutes the existing marks, or resamples
  from the existing marks.
}
\usage{
   rlabel(X, labels=marks(X), permute=TRUE)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"},
    \code{"lpp"}, \code{"pp3"} or \code{"ppx"}).
  }
  \item{labels}{
    Vector of values from which the new marks will be drawn
    at random. Defaults to the vector of existing marks.
  }
  \item{permute}{
    Logical value indicating whether to generate new marks
    by randomly permuting \code{labels} or
    by drawing a random sample with replacement.
  }
}
\value{
  A marked point pattern (of the same class as \code{X}).
}
\details{
  This very simple function allocates random marks to
  an existing point pattern \code{X}. It is useful
  for hypothesis testing purposes.

  In the simplest case, the command \code{rlabel(X)} yields
  a point pattern obtained from \code{X} by randomly permuting
  the marks of the points.

  If \code{permute=TRUE}, then \code{labels} should be a vector of
  length equal to the number of points in \code{X}.
  The result of \code{rlabel} will be a point pattern
  with locations given by \code{X} and marks given by
  a random permutation of \code{labels} (i.e. a random sample without
  replacement).

  If \code{permute=FALSE}, then \code{labels} may be a vector of
  any length. 
  The result of \code{rlabel} will be a point pattern
  with locations given by \code{X} and marks given by
  a random sample from \code{labels} (with replacement).
}
\seealso{
  \code{\link{marks<-}} to assign arbitrary marks.
}
\examples{
   data(amacrine)

   # Randomly permute the marks "on" and "off"
   # Result always has 142 "off" and 152 "on"
   Y <- rlabel(amacrine)

   # randomly allocate marks "on" and "off"
   # with probabilities p(off) = 0.48, p(on) = 0.52
   Y <- rlabel(amacrine, permute=FALSE)

   # randomly allocate marks "A" and "B" with equal probability
   data(cells)
   Y <- rlabel(cells, labels=factor(c("A", "B")), permute=FALSE)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
