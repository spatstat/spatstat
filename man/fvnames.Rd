\name{fvnames}
\alias{fvnames}
\alias{fvnames<-}
\title{
  Abbreviations for Groups of Columns in Function Value Table
}
\description{
  Groups of columns in a function value table (object
  of class \code{"fv"}) identified by standard abbreviations.
}
\usage{
fvnames(X, a = ".")

fvnames(X, a = ".") <- value

}
\arguments{
  \item{X}{
    Function value table (object of class \code{"fv"}).
    See \code{\link{fv.object}}.
  }
  \item{a}{
    One of the standard abbreviations listed below.
  }
  \item{value}{
    Character vector containing names of columns of \code{X}.
  }
}
\details{
  An object of class \code{"fv"} represents a table of
  values of a function, usually a summary function for spatial data
  such as the \eqn{K}-function, for which several different statistical
  estimators may be available. The different estimates are stored
  as columns of the table.

  Auxiliary information carried in the object \code{X} specifies some
  columns or groups of columns of this table that should be
  used for particular purposes.
  For convenience these groups can be referred to by standard
  abbreviations which are recognised by various functions
  in the \pkg{spatstat} package, such as \code{\link{plot.fv}}.

  These abbreviations are:
  \tabular{ll}{
    \code{".x"} \tab the function argument \cr
    \code{".y"} \tab the recommended value of the function \cr
    \code{"."}  \tab all function values to be plotted by default \cr
                \tab (in order of plotting) \cr
    \code{".s"} \tab the upper and lower limits of shading \cr
                \tab (for envelopes and confidence intervals)\cr
    \code{".a"} \tab all function values
  }
  The command \code{fvnames(X, a)} expands the abbreviation \code{a} and returns
  a character vector containing the names of the columns.

  The assignment \code{fvnames(X, a) <- value} changes the
  definition of the abbreviation \code{a} to the character vector
  \code{value}. It does not change the labels of any columns.

  Note that \code{fvnames(x, ".")} lists the columns of values that will
  be plotted by default, in the order that they would be plotted, not in
  order of the column position. The order in which curves are plotted
  affects the colours and line styles associated with the curves.
}
\value{
  For \code{fvnames}, a character vector.

  For \code{fvnames<-}, the updated object.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{fv.object}},
  \code{\link{plot.fv}}
}
\examples{
   K <- Kest(cells)
   fvnames(K, ".y")
   fvnames(K, ".y") <- "trans"
}
\keyword{spatial}
\keyword{manip}
