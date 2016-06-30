\name{is.multitype.ppm}
\alias{is.multitype.ppm}
\alias{is.multitype.lppm}
\title{Test Whether A Point Process Model is Multitype}
\description{
  Tests whether a fitted point process model involves ``marks''
  attached to the points that classify the points into several types.
}
\usage{
  \method{is.multitype}{ppm}(X, \dots) 

  \method{is.multitype}{lppm}(X, \dots) 
}
\arguments{
  \item{X}{
    Fitted point process model (object of class \code{"ppm"})
    usually obtained from \code{\link{ppm}}.
    Alternatively a model of class \code{"lppm"}.
  }
  \item{\dots}{
    Ignored.
  }
}
\value{
  Logical value, equal to \code{TRUE} if
  \code{X} is a model that was fitted to a multitype point pattern dataset.
}
\details{
  ``Marks'' are observations attached to each point of a point pattern.
  For example the \code{\link{longleaf}} dataset contains the locations
  of trees, each tree being marked by its diameter;
  the \code{\link{amacrine}} dataset gives the locations of cells
  of two types (on/off) and the type of cell may be regarded as a mark attached
  to the location of the cell.

  The argument \code{X} is a fitted point process model
  (an object of class \code{"ppm"}) typically obtained
  by fitting a model to point pattern data using \code{\link{ppm}}.

  This function returns \code{TRUE} if the \emph{original data}
  (to which the model \code{X} was fitted) were a multitype point pattern.

  Note that this is not the same as testing whether the
  model involves terms that depend on the marks (i.e. whether the
  fitted model ignores the marks in the data).
  Currently we have not implemented a test for this.

  If this function returns \code{TRUE}, the implications are
  (for example) that
  any simulation of this model will require simulation of random marks
  as well as random point locations.
}
\seealso{
  \code{\link{is.multitype}},
  \code{\link{is.multitype.ppp}}
}
\examples{
   X <- lansing
   # Multitype point pattern --- trees marked by species

   \testonly{
      # Smaller dataset
      X <- amacrine
   }

  fit1 <- ppm(X, ~ marks, Poisson())
  is.multitype(fit1)
  # TRUE

  fit2 <- ppm(X, ~ 1, Poisson())
  is.multitype(fit2)
  # TRUE

  # Unmarked point pattern
  fit3 <- ppm(cells, ~ 1, Poisson())
  is.multitype(fit3)
  # FALSE

}

\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{manip}
\keyword{models}
