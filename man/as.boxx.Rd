\name{as.boxx}
\alias{as.boxx}
\title{Convert Data to Multi-Dimensional Box}
\description{Interprets data as the dimensions of a multi-dimensional box.}
\usage{
  as.boxx(\dots, warn.owin = TRUE)
}
\arguments{
  \item{\dots}{
    Data that can be interpreted as giving the dimensions of a
    multi-dimensional box. See Details.
  }
  \item{warn.owin}{
    Logical value indicating whether to print a warning
    if a non-rectangular window (object of class \code{"owin"})
    is supplied.
  }
}
\details{
  Either a single argument should be provided which is one
  of the following:
  \itemize{
    \item an object of class \code{"boxx"}
    \item an object of class \code{"box3"}
    \item an object of class \code{"owin"}
    \item a numeric vector of even length, specifying the corners
    of the box. See Examples
  }
  or a list of arguments acceptable to \code{\link{boxx}}.
}
\value{A \code{"boxx"} object.}
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
 # Convert unit square to two dimensional box.
 W <- owin()
 as.boxx(W)
 # Make three dimensional box [0,1]x[0,1]x[0,1] from numeric vector
 as.boxx(c(0,1,0,1,0,1))
}
