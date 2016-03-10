\name{Kmodel}
\alias{Kmodel}
\alias{pcfmodel}
\title{K Function or Pair Correlation Function of a Point Process Model}
\description{
  Returns the theoretical \eqn{K} function or the pair correlation function
  of a point process model.
}
\usage{
   Kmodel(model, \dots)

   pcfmodel(model, \dots)
}
\arguments{
  \item{model}{
    A fitted point process model of some kind.
  }
  \item{\dots}{
    Ignored.
  }
}
\value{
  A \code{function} in the \R language,
  which takes one argument \code{r}.
}
\details{
  For certain types of point process models, it is possible to
  write down a mathematical expression for the \eqn{K} function
  or the pair correlation function of the model. 
  
  The functions \code{Kmodel} and \code{pcfmodel} give the
  theoretical \eqn{K}-function and the theoretical pair correlation
  function for a point process model that has been fitted to data.
  
  The functions \code{Kmodel} and \code{pcfmodel} are generic,
  with methods for the classes \code{"kppm"} (cluster processes
  and Cox processes) and \code{"ppm"} (Gibbs processes).
  
  The return value is a \code{function} in the \R language,
  which takes one argument \code{r}.
  Evaluation of this function, on a numeric vector \code{r},
  yields values of the desired \eqn{K} function or pair correlation
  function at these distance values.
}
\seealso{
  \code{\link{Kest}} or \code{\link{pcf}}
  to estimate the \eqn{K} function or pair correlation function
  nonparametrically from data.

  \code{\link{Kmodel.kppm}} for the method for cluster processes
  and Cox processes.
  
  \code{\link{Kmodel.ppm}} for the method for Gibbs processes.
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{models}
