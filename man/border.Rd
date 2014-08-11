\name{border}
\alias{border}
\title{Border Region of a Window}
\description{
  Computes the border region of a window,
  that is, the region lying within a specified distance of the boundary
  of a window.
}
\usage{
border(w, r, outside=FALSE, ...)
}
\arguments{
  \item{w}{A window (object of class \code{"owin"})
    or something acceptable to \code{\link{as.owin}}.
  }
  \item{r}{Numerical value.}
  \item{outside}{Logical value determining whether to compute the
    border outside or inside \code{w}.}
  \item{\dots}{
    Optional arguments passed to \code{\link{erosion}}
    (if \code{outside=FALSE}) or to \code{\link{dilation}}
    (if \code{outside=TRUE}).
  }
}
\value{
  A window (object of class \code{"owin"}).
}
\details{
  By default (if \code{outside=FALSE}),
  the border region is the subset of \code{w}
  lying within a distance \code{r} of the boundary of \code{w}.
  It is computed by eroding \code{w} by the distance \code{r} (using
  \code{\link{erosion}}) and
  subtracting this eroded window from the original window \code{w}.

  If \code{outside=TRUE}, the border region is the set of locations
  outside \code{w} lying within a distance \code{r} of \code{w}.
  It is computed by dilating \code{w} by the distance \code{r}
  (using \code{\link{dilation}}) and
  subtracting the original window \code{w} from the dilated window.
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{erosion}},
  \code{\link{dilation}}
}
\examples{
# rectangle
   u <- unit.square()
   border(u, 0.1)
   border(u, 0.1, outside=TRUE)
# polygon
   \testonly{opa <- spatstat.options(npixel=32)}
   data(letterR)
   plot(letterR)
   plot(border(letterR, 0.1), add=TRUE)
   plot(border(letterR, 0.1, outside=TRUE), add=TRUE)
   \testonly{spatstat.options(opa)}
}
\keyword{spatial}
\keyword{math}
