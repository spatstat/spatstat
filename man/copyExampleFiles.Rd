\name{copyExampleFiles}
\alias{copyExampleFiles}
\title{
  Copy Data Files for Example 
}
\description{
  This command copies several data files to
  a folder (directory) chosen by the user,
  so that they can be used for a practice example.
}
\usage{
copyExampleFiles(which, folder = getwd())
}
\arguments{
  \item{which}{
    Character string name (partially matched)
    of one of the datasets installed in \code{spatstat} for which
    the original data files are provided. 
    If \code{which} is missing, a list of available options is printed.
  }
  \item{folder}{
    Character string path name of a folder (directory) in which the
    files will be placed. Defaults to the current working directory.
  }
}
\details{
  The original text files containing data for the selected
  dataset are copied to the chosen folder.

  This is part of an exercise described in Chapter 3 of
  Baddeley, Rubak and Turner (2015).
}
\references{
  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Analysing Spatial Point Patterns with R}.
  Chapman and Hall/CRC, to appear.
}
\author{
\adrian


\rolf

and \ege

}
\examples{
   copyExampleFiles()
}
\seealso{
  \code{\link{vesicles}},
  \code{\link{finpines}}
}
\keyword{spatial}
