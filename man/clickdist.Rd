\name{clickdist}
\alias{clickdist}
\title{Interactively Measure Distance}
\description{
  Measures the distance between two points
  which the user has clicked on.
}
\usage{
  clickdist()
}
\value{
  A single nonnegative number.
}
\details{
  This function allows the user to measure the distance
  between two spatial locations, interactively,
  by clicking on the screen display.

  When \code{clickdist()} is called, the user is expected to
  click two points in the current graphics device. The distance
  between these points will be returned.

  This function uses the \R{} command \code{\link[graphics]{locator}} to
  input the mouse clicks. It only works on screen devices such as
  \sQuote{X11}, \sQuote{windows} and \sQuote{quartz}.
}
\seealso{
  \code{\link[graphics]{locator}},
  \code{\link{clickppp}},
  \code{\link{clickpoly}},
  \code{\link{clickbox}}
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\keyword{spatial}
\keyword{iplot}
