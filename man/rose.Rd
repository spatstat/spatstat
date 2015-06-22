\name{rose}
\alias{rose}
\alias{rose.default}
\alias{rose.histogram}
\alias{rose.density}
\alias{rose.fv}
\title{Rose Diagram}
\description{
  Plots a rose diagram (rose of directions),
  the analogue of a histogram or density plot for angular data.
}
\usage{
rose(x, \dots)

\method{rose}{default}(x, breaks = NULL, \dots,
                       weights=NULL,
                       nclass = NULL,
                       unit = c("degree", "radian"),
                       start=0, clockwise=FALSE,
                       main)

\method{rose}{histogram}(x, \dots,
                       unit = c("degree", "radian"),
                       start=0, clockwise=FALSE,
                       main, labels=TRUE, at=NULL, do.plot = TRUE)

\method{rose}{density}(x, \dots, 
                  unit = c("degree", "radian"), 
                  start=0, clockwise=FALSE,
                  main, labels=TRUE, at=NULL, do.plot = TRUE)

\method{rose}{fv}(x, \dots, 
                  unit = c("degree", "radian"),
                  start=0, clockwise=FALSE,
                  main, labels=TRUE, at=NULL, do.plot = TRUE)
}
\arguments{
  \item{x}{
    Data to be plotted.
    A numeric vector containing angles,
    or a \code{histogram} object containing a histogram of
    angular values, or a \code{density} object containing a
    smooth density estimate for angular data,
    or an \code{fv} object giving a function
    of an angular argument.
  }
  \item{breaks, nclass}{
    Arguments passed to \code{\link[graphics]{hist}} to determine
    the histogram breakpoints.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link[graphics]{polygon}}
    controlling the appearance of the plot
    (or passed from \code{rose.default} to \code{\link[graphics]{hist}}
    to control the calculation of the histogram).
  }
  \item{unit}{
    The unit in which the angles are expressed.
  }
  \item{start}{
    The starting direction for measurement of angles,
    that is, the spatial direction which corresponds to a measured angle
    of zero. Either a character string giving a compass direction
    (\code{"N"} for north, \code{"S"} for south, \code{"E"} for east,
    or \code{"W"} for west) or a number giving the angle from the
    the horizontal (East) axis to the starting direction. For example,
    if \code{unit="degree"} and \code{clockwise=FALSE},
    then \code{start=90} and \code{start="N"} are equivalent.
    The default is to measure angles anti-clockwise from the
    horizontal axis (East direction).
  }
  \item{clockwise}{
    Logical value indicating whether angles increase in the clockwise
    direction (\code{clockwise=TRUE}) or anti-clockwise,
    counter-clockwise direction (\code{clockwise=FALSE}, the default).
  }
  \item{weights}{
    Optional vector of numeric weights associated with \code{x}.
  }
  \item{main}{
    Optional main title for the plot.
  }
  \item{labels}{
    Either a logical value indicating whether to plot labels
    next to the tick marks, or a vector of labels for the tick marks.
  }
  \item{at}{
    Optional vector of angles at which tick marks should be plotted.
    Set \code{at=numeric(0)} to suppress tick marks.
  }
  \item{do.plot}{
    Logical value indicating whether to really perform the plot.
  }
}
\details{
  A rose diagram or rose of directions is the analogue of
  a histogram or bar chart for data which represent angles
  in two dimensions. The bars of the bar chart are replaced by
  circular sectors in the rose diagram.

  The function \code{rose} is generic, with a default method
  for numeric data, and methods for histograms and function tables.
  
  If \code{x} is a numeric vector, it must contain angular values
  in the range 0 to 360 (if \code{unit="degree"})
  or in the range 0 to \code{2 * pi} (if \code{unit="radian"}).
  A histogram of the data will first be computed using
  \code{\link[graphics]{hist}}. Then the rose diagram of this histogram
  will be plotted by \code{rose.histogram}.

  If \code{x} is an object of class \code{"histogram"} produced by
  the function \code{\link[graphics]{hist}}, representing the histogram
  of angular data, then the rose diagram of the densities
  (rather than the counts) in this histogram object will be plotted. 

  If \code{x} is an object of class \code{"density"} produced by
  \code{\link{circdensity}} or \code{\link[stats]{density.default}},
  representing a kernel smoothed density estimate of angular data,
  then the rose diagram of the density estimate will be plotted. 
  
  If \code{x} is a function value table (object of class \code{"fv"})
  then the argument of the function will be interpreted as an angle,
  and the value of the function will be interpreted as the radius.

  By default, angles are interpreted using the mathematical convention
  where the zero angle is the horizontal \eqn{x} axis, and angles
  increase anti-clockwise. Other conventions can be specified
  using the arguments \code{start} and \code{clockwise}.
  Standard compass directions are obtained by setting \code{unit="degree"},
  \code{start="N"} and \code{clockwise=TRUE}.  
}
\value{A window (class \code{"owin"}) containing the plotted region.}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\seealso{
  \code{\link{fv}}, \code{\link[graphics]{hist}},
  \code{\link{circdensity}},
  \code{\link[stats]{density.default}}.
}
\examples{
  ang <- runif(1000, max=360)
  rose(ang, col="grey")
  rose(ang, col="grey", start="N", clockwise=TRUE)
}
\keyword{spatial}
\keyword{hplot}

