\name{methods.ppx}
\Rdversion{1.1}
\alias{methods.ppx} %DoNotExport
\alias{print.ppx}
\alias{plot.ppx}
\alias{unitname.ppx}
\alias{unitname<-.ppx}
\alias{scale.ppx}
\title{
  Methods for Multidimensional Space-Time Point Patterns
}
\description{
  Methods for printing and plotting a general multidimensional
  space-time point pattern.
}
\usage{
\method{print}{ppx}(x, ...)
\method{plot}{ppx}(x, ...)
\method{unitname}{ppx}(x)
\method{unitname}{ppx}(x) <- value
\method{scale}{ppx}(x, center=TRUE, scale=TRUE)
}
\arguments{
  \item{x}{
    Multidimensional point pattern (object of class \code{"ppx"}).
  }
  \item{\dots}{
    Additional arguments passed to plot methods.
  }
  \item{value}{
    Name of the unit of length. See \code{\link{unitname}}.
  }
  \item{center,scale}{
    Arguments passed to \code{\link[base]{scale.default}}
    to determine the rescaling.
  }
}
\details{
  These are methods for the generic functions
  \code{\link{print}}, \code{\link{plot}},
  \code{\link{unitname}},  \code{\link{unitname<-}}
  and \code{\link[base]{scale}}
  for the class \code{"ppx"} of multidimensional point patterns.
  
  The \code{print} method prints a description
  of the point pattern and its spatial domain.

  The \code{unitname} method extracts the name of the unit of
  length in which the point coordinates are expressed.
  The \code{unitname<-} method assigns the name of the unit of length.

  The \code{scale} method rescales each spatial coordinate of \code{x}.
}
\value{
  For \code{print.ppx} and \code{plot.ppx} the value is \code{NULL}.
  For \code{unitname.ppx} the value is an object of class \code{"units"}.
  For \code{unitname<-.ppx} and \code{scale.ppx}
  the value is another object of class \code{"ppx"}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{ppx}},
  \code{\link{unitname}}  
}
\keyword{spatial}
