\name{collapse.fv}
\alias{collapse.fv}
\alias{collapse.anylist}
\title{
  Collapse Several Function Tables into One
}
\description{
  Combines several function tables (objects of class \code{"fv"})
  into a single function table, merging columns that are identical
  and relabelling columns that are different.
}
\usage{
\method{collapse}{fv}(object, \dots, same = NULL, different = NULL)

\method{collapse}{anylist}(object, \dots, same = NULL, different = NULL)
}
\arguments{
  \item{object}{
    An object of class \code{"fv"}, or a list of such objects.
  }
  \item{\dots}{
    Additional objects of class \code{"fv"}.
  }
  \item{same}{
    Character string or character vector specifying a column or columns,
    present in each \code{"fv"} object, that are identical
    in each object. This column or columns will be included only once.
  }
  \item{different}{
    Character string or character vector specifying a column or columns,
    present in each \code{"fv"} object, that contain different values in
    each object. Each of these columns of data will be included, with
    labels that distinguish them from each other.
  }
}
\details{
  This is a method for the generic function
  \code{\link[nlme]{collapse}}.
  
  It combines the data in several function tables 
  (objects of class \code{"fv"}, see \code{\link{fv.object}})
  to make a single function table.
  It is essentially a smart wrapper for
  \code{\link{cbind.fv}}.

  A typical application is to calculate the same summary statistic
  (such as the \eqn{K} function) for different point patterns,
  and then to use \code{collapse.fv} to combine the results into a
  single object that can easily be plotted. See the Examples.
  
  The arguments \code{object} and \code{\dots} should be function tables
  (objects of class \code{"fv"}, see \code{\link{fv.object}})
  that are compatible in the sense that they
  have the same values of the function argument. 

  The argument \code{same} identifies any columns that are present
  in each function table, and which are known to contain exactly
  the same values in each table. This column or columns will be
  included only once in the result.

  The argument \code{different} identifies any columns that are present
  in each function table, and which contain different numerical values
  in each table. Each of these columns will be included, with labels
  to distinguish them. 

  Columns that are not named in \code{same} or \code{different} will not
  be included.
}
\value{
  Object of class \code{"fv"}.
}
\seealso{
  \code{\link{fv.object}},
  \code{\link{cbind.fv}}
}
\examples{
  # generate simulated data
  X <- replicate(3, rpoispp(100), simplify=FALSE)
  names(X) <- paste("Simulation", 1:3)
  # compute K function estimates
  Klist <- anylapply(X, Kest)
  # collapse
  K <- collapse(Klist, same="theo", different="iso")
  K
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
