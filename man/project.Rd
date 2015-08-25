\name{project}
\alias{project}
\title{
  Force Model to be Valid
}
\description{
  Check whether a model is valid, and if not,
  find the nearest model which is valid.
}
\usage{
project(object, \dots)
}
\arguments{
  \item{object}{
    A statistical model, belonging to some class.
  }
  \item{\dots}{Arguments passed to methods.}
}
\details{
  The function \code{project} is generic,
  and has methods for several classes of statistical models
  in the \pkg{spatstat} package (mostly point process models).
  Its purpose is to check whether a given model is valid
  (for example, that none of the model parameters are \code{NA})
  and, if not, to find the nearest model which is valid.

  See the methods for more information.
}
\value{
  Another model of the same kind.
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
  \code{\link{project.ppm}},
  \code{\link{project.lppm}},
  \code{\link{valid}}.
}
\keyword{spatial}
\keyword{models}
