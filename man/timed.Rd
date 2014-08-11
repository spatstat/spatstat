\name{timed}
\alias{timed}
\title{
  Record the Computation Time
}
\description{
  Saves the result of a calculation as an object of class \code{"timed"}
  which includes information about the time taken to compute the result.
  The computation time is printed when the object is printed.
}
\usage{
timed(x, ..., starttime = NULL, timetaken = NULL)
}
\arguments{
  \item{x}{
    An expression to be evaluated, or an object that has already
    been evaluated.
  }
  \item{starttime}{
    The time at which the computation is defined to have started.
    The default is the current time.
    Ignored if \code{timetaken} is given.
  }
  \item{timetaken}{
    The length of time taken to perform the computation.
    The default is the time taken to evaluate \code{x}.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a simple mechanism for recording how long it takes to
  perform complicated calculations (usually for the purposes of
  reporting in a publication).
  
  If \code{x} is an expression to be evaluated, \code{timed(x)}
  evaluates the expression and measures the
  time taken to evaluate it.
  The result is saved as an object of the class
  \code{"timed"}. Printing this object displays the computation time.

  If \code{x} is an object which has already been computed,
  then the time taken to compute the object can be specified either directly
  by the argument \code{timetaken}, or indirectly by the argument
  \code{starttime}.
  \itemize{
    \item
    \code{timetaken} is the duration of time taken to perform
    the computation. It should be the difference of two clock times
    returned by \code{\link{proc.time}}. Typically the user
    sets \code{begin <- proc.time()} before commencing the calculations,
    then \code{end <- proc.time()} after completing the calculations,
    and then sets \code{timetaken <- end - begin}.
    \item
    \code{starttime} is the clock time at which the computation started.
    It should be a value that was returned by \code{\link{proc.time}}
    at some earlier time when the calculations commenced.
    When \code{timed} is called, the computation time will be taken
    as the difference between the current clock time and
    \code{starttime}. Typically the user
    sets \code{begin <- proc.time()} before commencing the calculations,
    and when the calculations are completed, the user calls
    \code{result <- timed(result, starttime=begin)}. 
  }
  
  If the result of evaluating \code{x} belongs to other S3 classes,
  then the result of \code{timed(x, \dots)} also inherits these classes,
  and printing the object will display the appropriate information for these
  classes as well.
}
\value{
  An object inheriting the class \code{"timed"}.
}
\examples{
  timed(envelope(cells, nsim=19))  
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
