\name{vesiclesCopy}
\alias{vesiclesCopy}
\title{
  Copy the Data Files for Synaptic Vesicles Example 
}
\description{
  This command copies several data files to
  a folder (directory) chosen by the user,
  so that they can be used for a practice example.
}
\usage{
vesiclesCopy(folder = getwd())
}
\arguments{
  \item{folder}{
    Character string name of a folder (directory) in which the
    files will be placed. Defaults to the current working directory.
  }
}
\details{
  The original text files containing data for the
  \code{\link{vesicles}} dataset are copied to
  the chosen folder.

  This is part of an exercise described in Chapter 3 of
  Baddeley, Rubak and Turner (2015).
}
\references{
  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Analysing Spatial Point Patterns with R}.
  Chapman and Hall/CRC, to appear.
}
\author{
Adrian Baddeley
\email{Adrian.Baddeley@uwa.edu.au}
\url{http://www.maths.uwa.edu.au/~adrian/}
Rolf Turner
\email{r.turner@auckland.ac.nz}
and Ege Rubak
\email{rubak@math.aau.dk}
}
\seealso{
  \code{\link{vesicles}}
}
\keyword{spatial}
