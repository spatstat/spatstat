\name{print.ppp}
\alias{print.ppp}
\title{Print Brief Details of a Point Pattern Dataset}
\description{
  Prints a very brief description of a point pattern dataset.
}
\usage{
  \method{print}{ppp}(x, \dots)
}
\arguments{
  \item{x}{Point pattern (object of class \code{"ppp"}).}
  \item{\dots}{Ignored.}
}
\details{
  A very brief description of the point pattern \code{x} is printed.

  This is a method for the generic function \code{\link{print}}.
}
\seealso{
  \code{\link{print}},
  \code{\link{print.owin}},
  \code{\link{summary.ppp}}
}
\examples{
  data(cells)      # plain vanilla point pattern
  cells

  data(lansing)    # multitype point pattern
  lansing          
  
  data(longleaf)    # numeric marks
  longleaf          

  data(demopat)     # weird polygonal window
  demopat
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{print}

