\name{ppp.object}
\alias{ppp.object} %DoNotExport
\title{Class of Point Patterns}
\description{
  A class \code{"ppp"} to represent a two-dimensional point
  pattern. Includes information about the window in which the
  pattern was observed. Optionally includes marks.
}
\details{
  This class represents
  a two-dimensional point pattern dataset. It specifies
  \itemize{
    \item the locations of the points
    \item the window in which the pattern was observed
    \item optionally, ``marks'' attached to each point
      (extra information such as a type label).
    }
  If \code{X} is an object of type \code{ppp},
  it contains the following elements:
  \tabular{ll}{
    \code{x} \tab vector of \eqn{x} coordinates of data points \cr
    \code{y} \tab vector of \eqn{y} coordinates of data points \cr
    \code{n} \tab number of points \cr
    \code{window} \tab window of observation \cr
                  \tab (an object of class \code{\link{owin}}) \cr
    \code{marks} \tab optional vector or data frame of marks
  }
  Users are strongly advised not to manipulate these entries
  directly.
  
  Objects of class \code{"ppp"}
  may be created by the function
  \code{\link{ppp}}
  and converted from other types of data by the function
  \code{\link{as.ppp}}.
  Note that you must always specify the window of observation;
  there is intentionally no default action of ``guessing'' the window
  dimensions from the data points alone.

  Standard point pattern datasets provided with the package
  include
  \code{\link{amacrine}},
  \code{\link{betacells}},
  \code{\link{bramblecanes}},
  \code{\link{cells}},
  \code{\link{demopat}},
  \code{\link{ganglia}},
  \code{\link{lansing}},
  \code{\link{longleaf}},
  \code{\link{nztrees}},
  \code{\link{redwood}},
  \code{\link{simdat}} and
  \code{\link{swedishpines}}.
  
  Point patterns may be scanned from your own data files by
  \code{\link{scanpp}} or by using \code{\link{read.table}} and
  \code{\link{as.ppp}}.
  
  They may be manipulated by the functions 
  \code{\link{[.ppp}} and
  \code{\link{superimpose}}.

  Point pattern objects can be plotted just by typing \code{plot(X)}
  which invokes the \code{plot} method for point pattern objects,
  \code{\link{plot.ppp}}. See \code{\link{plot.ppp}} for further information.

  There are also methods for \code{summary} and \code{print}
  for point patterns. Use \code{summary(X)} to see a useful description
  of the data.
  
  Patterns may be generated at random by
   \code{\link{runifpoint}},
   \code{\link{rpoispp}},
   \code{\link{rMaternI}},
   \code{\link{rMaternII}},
   \code{\link{rSSI}},
   \code{\link{rNeymanScott}},
   \code{\link{rMatClust}},
   and
   \code{\link{rThomas}}.

   Most functions which are intended to operate on a window
   (of class \code{\link{owin}})
   will, if presented with a \code{\link{ppp}} object instead,
   automatically extract the window information from the point pattern.
}
\seealso{
  \code{\link{owin}},
  \code{\link{ppp}},
  \code{\link{as.ppp}},
  \code{\link{[.ppp}}
}
\section{Warnings}{
  The internal representation of marks is likely to change in the
  next release of this package.
}
\examples{
  x <- runif(100)
  y <- runif(100)
  X <- ppp(x, y, c(0,1),c(0,1))
  X
  \dontrun{plot(X)}
  mar <- sample(1:3, 100, replace=TRUE)
  mm <- ppp(x, y, c(0,1), c(0,1), marks=mar)
  \dontrun{plot(mm)}
  # points with mark equal to 2
  ss <- mm[ mm$marks == 2 , ]
  \dontrun{plot(ss)}
  # left half of pattern 'mm'
  lu <- owin(c(0,0.5),c(0,1))
  mmleft <- mm[ , lu]
  \dontrun{plot(mmleft)}
  \dontrun{
  # input data from file
  qq <- scanpp("my.table", unit.square())

  # interactively build a point pattern
  plot(unit.square())
  X <- as.ppp(locator(10), unit.square())
  plot(X)
  }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{attribute}
 
 
