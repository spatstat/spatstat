\name{progressreport}
\alias{progressreport}
\title{Print Progress Reports}
\description{
  Prints Progress Reports during a loop or iterative calculation.
}
\usage{
progressreport(i, n,
               every = min(100,max(1, ceiling(n/100))),
               nperline = min(charsperline,
                              every * ceiling(charsperline/(every + 3))),
               charsperline = 60,
               style=spatstat.options("progress"),
               state=NULL)
}
\arguments{
  \item{i}{
    Integer. The current iteration number (from 1 to \code{n}).
  }
  \item{n}{
    Integer. The (maximum) number of iterations to be computed.
  }
  \item{every}{
    Optional integer. The number of iterations between successive
    reports.
  }
  \item{nperline}{
    Optional integer. The maximum number of reports to be printed per
    line of output.
  }
  \item{charsperline}{
    Optional integer. The number of characters in a line of output.
  }
  \item{style}{
    Character string determining the style of display.
    See Details.
  }
  \item{state}{
    Optional. A list containing the internal data.
  }
}
\details{
  This is a convenient function for reporting progress
  during an iterative sequence of calculations
  or a suite of simulations.

  If \code{style="tk"} then \code{tcltk::tkProgressBar} is
  used to pop-up a new graphics window showing a progress bar.
  This requires the package \pkg{tcltk}.

  If \code{style="txtbar"} then \code{\link[utils]{txtProgressBar}} is
  used to represent progress as a bar made of text characters in the
  \R interpreter window.

  If \code{style="tty"}, then progress reports are printed using
  \code{cat}. This only seems to work well under Linux.

  It is optional, but strongly advisable, to use the argument \code{state}
  to store and update the internal data for the progress reports
  (such as the cumulative time taken for computation)
  as shown in the last example below.
  This avoids conflicts with other programs that might be
  calling \code{progressreport} at the same time.
}
\value{
  If \code{state} was \code{NULL}, the result is \code{NULL}.
  Otherwise the result is the updated value of \code{state}.
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
  for(i in 1:40) {
     #
     # code that does something...
     # 
     progressreport(i, 40)
  }

  # saving internal state
  sta <- list()
  for(i in 1:20) {
     # some code ...
     sta <- progressreport(i, 20, state=sta)
  }
}
\keyword{print}
