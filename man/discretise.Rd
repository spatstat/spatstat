\name{discretise}
\alias{discretise}
\title{
  Safely Convert Point Pattern Window to Binary Mask
}
\description{
  Given a point pattern, discretise its window by converting it to a
  binary pixel mask, adjusting the mask so that it still contains all
  the points.
}
\usage{
  discretise(X, eps = NULL, dimyx = NULL, xy = NULL)
}
\arguments{
  \item{X}{A point pattern (object of class \code{"ppp"}) to be converted.}
  \item{eps}{(optional) width and height of each pixel}
  \item{dimyx}{(optional) pixel array dimensions}
  \item{xy}{(optional) pixel coordinates}
}
\details{
  This function modifies the point pattern \code{X} by converting its
  observation window \code{Window(X)} to a binary pixel image (a window
  of type \code{"mask"}). It ensures that no points of \code{X} are
  deleted by the discretisation.

  The window is first discretised using \code{\link{as.mask}}. 
  It can happen that points of \code{X} that were inside the original
  window may fall outside the new mask.
  The \code{discretise} function corrects this by augmenting the mask
  (so that the mask includes any pixel that contains a point of the pattern).

  The arguments \code{eps}, \code{dimyx} and \code{xy}
  control the fineness of the pixel array. They are passed to
  \code{\link{as.mask}}.
  
  If \code{eps}, \code{dimyx} and \code{xy} are all absent or
  \code{NULL}, and if the window of \code{X} is of type \code{"mask"}
  to start with, then \code{discretise(X)} returns \code{X} unchanged.

  See \code{\link{as.mask}} for further details
  about the arguments \code{eps}, \code{dimyx},
  and \code{xy}, and the process of converting
  a window to one of type \code{mask}.
}
\section{Error checking}{
  Before doing anything, \code{discretise} checks that
  all the points of the pattern are actually
  inside the original window.  This is guaranteed to
  be the case if the pattern was constructed using \code{\link{ppp}}
  or \code{\link{as.ppp}}. However anomalies are possible if the
  point pattern was created or manipulated inappropriately.
  These will cause an error.
}
\value{
  A point pattern (object of class \code{"ppp"}),
  identical to \code{X}, except that
  its observation window has been converted to one
  of type \code{mask}. 
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\seealso{
  \code{\link{as.mask}}
}
\examples{
  data(demopat)
  X <- demopat
  plot(X, main="original pattern")
  Y <- discretise(X, dimyx=50)
  plot(Y, main="discretise(X)")
  stopifnot(npoints(X) == npoints(Y))

  # what happens if we just convert the window to a mask?
  W <- Window(X)
  M <- as.mask(W, dimyx=50)
  plot(M, main="window of X converted to mask")
  plot(X, add=TRUE, pch=16)
  plot(X[M], add=TRUE, pch=1, cex=1.5)
  XM <- X[M]
  cat(paste(npoints(X) - npoints(XM), "points of X lie outside M\n"))
}
\keyword{spatial}
\keyword{manip}


