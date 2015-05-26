\name{begins}
\alias{begins}
\title{
  Check Start of Character String 
}
\description{
  Checks whether a character string begins with a particular prefix.
}
\usage{
begins(x, firstbit)
}
\arguments{
  \item{x}{
    Character string, or vector of character strings, to be tested.
  }
  \item{firstbit}{
    A single character string.
  }
}
\details{
  This simple wrapper function checks whether (each entry in) \code{x}
  begins with the string \code{firstbit}, and returns a logical value
  or logical vector with one entry for each entry of \code{x}.
  This function is useful mainly for reducing complexity in model formulae.
}
\value{
  Logical vector of the same length as \code{x}.
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
\examples{
  begins(c("Hello", "Goodbye"), "Hell")
  begins("anything", "")
}
\keyword{character}
