\name{timeTaken}
\alias{timeTaken}
\title{
  Extract the Total Computation Time
}
\description{
  Given an object or objects that contain timing information
  (reporting the amount of computer time taken to compute each object),
  this function extracts the timing data and evaluates the total time taken.
}
\usage{
timeTaken(..., warn=TRUE)
}
\arguments{
  \item{\dots}{
    One or more objects of class \code{"timed"} containing
    timing data. 
  }
  \item{warn}{
    Logical value indicating whether a warning should be issued if
    some of the arguments do not contain timing information.
  }
}
\details{
  An object of class \code{"timed"} contains information on
  the amount of computer time that was taken to compute the object.
  See \code{\link{timed}}.

  This function extracts the timing information from one or more
  such objects, and calculates the total time.
}
\value{
  An object inheriting the class \code{"timed"}.
}
\examples{
  A <- timed(Kest(cells))
  B <- timed(Gest(cells))
  A
  B
  timeTaken(A,B)
}
\seealso{
  \code{\link{timed}}
}
\author{
  \spatstatAuthors.
}
\keyword{utilities}
