\name{square}
\alias{square}
\alias{unit.square}
\title{Square Window}
\description{
Creates a square window
}
\usage{
 square(r=1, unitname=NULL)
 unit.square()
}
\arguments{
  \item{r}{Numeric. The side length of the square,
    or a vector giving the minimum and maximum coordinate values.
  }
  \item{unitname}{
    Optional. Name of unit of length. Either a single character string,
    or a vector of two character strings giving the
    singular and plural forms, respectively.
  }
}
\value{
  An object of class \code{"owin"} (see \code{\link{owin.object}})
  specifying a window. 
}
\details{
  If \code{r} is a number, \code{square(r)}
  is a shortcut for creating a window object
  representing the square 
  \eqn{[0,r] \times [0,r]}{[0,r] * [0,r]}.
  It is equivalent to the command
  \code{owin(c(0,r),c(0,r))}.

  If \code{r} is a vector of length 2, then
  \code{square(r)} creates the square with \code{x} and \code{y} coordinates
  ranging from \code{r[1]} to \code{r[2]}.
  
  \code{unit.square} creates the unit square
  \eqn{[0,1] \times [0,1]}{[0,1] * [0,1]}.
  It is equivalent to
  \code{square(1)} or \code{square()} or \code{owin(c(0,1),c(0,1))}.

  These commands are included for convenience,
  and to improve the readability of some code.
}
\seealso{
  \code{\link{owin.object}},
  \code{\link{owin}}
}
\examples{
 W <- square(10)
 W <- square(c(-1,1))
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{datagen}
 
 
