\name{is.subset.owin}
\alias{is.subset.owin}
\title{Determine Whether One Window is Contained In Another}
\description{
  Tests whether window \code{A} is a subset of window \code{B}.
}
\usage{
is.subset.owin(A, B)
}
\arguments{
  \item{A}{A window object (see Details).}
  \item{B}{A window object (see Details).}
}
\value{
  Logical scalar; \code{TRUE} if \code{A} is a sub-window of \code{B},
  otherwise \code{FALSE}.
}
\details{
  This function tests whether the window \code{A} is a subset
  of the window \code{B}.

  The arguments \code{A} and \code{B} must be window objects
  (either objects of class \code{"owin"}, or data that can be
  coerced to this class by \code{\link{as.owin}}).

  Various algorithms are used, depending on the geometrical type
  of the two windows.

  Note that if \code{B} is not rectangular, the algorithm proceeds by
  discretising \code{A}, converting it to a pixel mask using 
  \code{\link{as.mask}}. In this case the resulting
  answer is only ``approximately correct''. The accuracy of the
  approximation can be controlled: see \code{\link{as.mask}}.
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}

\examples{
   w1 <- as.owin(c(0,1,0,1))
   w2 <- as.owin(c(-1,2,-1,2))
   is.subset.owin(w1,w2)  # Returns TRUE.
   is.subset.owin(w2,w1)  # Returns FALSE.
}
\keyword{spatial}
\keyword{math}
