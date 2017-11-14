\name{reload.or.compute}
\alias{reload.or.compute}
\title{
  Compute Unless Previously Saved
}
\description{
  If the designated file does not yet exist,
  evaluate the expression and save the results in the file.
  If the file already exists, re-load the results from the file.
}
\usage{
reload.or.compute(filename, expr, objects = NULL,
                  destination = parent.frame(), force=FALSE)
}
\arguments{
  \item{filename}{
    Name of data file. A character string.
  }
  \item{expr}{
    \R language expression to be evaluated.
  }
  \item{objects}{
    Optional character vector of names of objects to be saved
    in \code{filename} after evaluating \code{expr},
    or names of objects that should be present in \code{filename}
    when loaded.
  }
  \item{destination}{
    Environment in which the resulting objects should be assigned.
  }
  \item{force}{
    Logical value indicating whether to perform the computation
    in any case.
  }
}
\details{
  This facility is useful for saving, and later re-loading, the results of
  time-consuming computations. It would typically be
  used in an \R script file or an \code{\link[utils]{Sweave}} document.

  If the file called \code{filename} does not yet exist,
  then \code{expr} will be evaluated
  and the results will be saved in \code{filename}.
  The optional argument \code{objects} specifies which results should be saved
  to the file: the default is to save all objects that were created
  by evaluating the expression. 

  If the file called \code{filename} already exists, then it
  will be loaded.
  The optional argument \code{objects} specifies the names
  of objects that should be present in the file; a warning is issued
  if any of them are missing.

  The resulting objects can be assigned into any desired \code{destination}.
  The default behaviour is equivalent to evaluating \code{expr}
  in the current environment.

  If \code{force=TRUE} then \code{expr} will be evaluated
  (regardless of whether the file already exists or not)
  and the results will be saved in \code{filename}, overwriting
  any previously-existing file with that name. This is a convenient
  way to force the code to re-compute everything
  in an \R script file or \code{\link[utils]{Sweave}} document.
}
\value{
  Character vector (invisible) giving the names of the objects computed
  or loaded.
}
\examples{
   \dontrun{
    if(FALSE) {
     reload.or.compute("mydata.rda", {
        x <- very.long.computation()
        y <- 42
      })
   }
   }
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{utilities}
