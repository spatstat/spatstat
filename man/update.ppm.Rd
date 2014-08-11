\name{update.ppm}
\alias{update.ppm}
\title{Update a Fitted Point Process Model}
\description{
  \code{update} method for class \code{"ppm"}.
}
\usage{
  \method{update}{ppm}(object, \dots, fixdummy=TRUE, use.internal=NULL,
                                      envir=environment(terms(object)))

}
\arguments{
  \item{object}{
    An existing fitted point process model,
    typically produced by \code{\link{ppm}}.
  }
  \item{\dots}{
    Arguments to be updated in the new call to \code{\link{ppm}}.
  }
  \item{fixdummy}{
    Logical flag indicating whether the quadrature scheme for the
    call to \code{\link{ppm}} should use the same set of dummy points
    as that in the original call.
  }
  \item{use.internal}{
    Optional. Logical flag indicating whether the model should be
    refitted using the internally saved data (\code{use.internal=TRUE})
    or by re-evaluating these data in the
    current frame (\code{use.internal=FALSE}).
  }
  \item{envir}{
    Environment in which to re-evaluate the call to \code{\link{ppm}}.
  }
}
\details{
  This is a method for the generic function \code{\link{update}}
  for the class \code{"ppm"}. An object of class \code{"ppm"}
  describes a fitted point process model. See \code{\link{ppm.object}})
  for details of this class.

  \code{update.ppm} will modify the point process model
  specified by \code{object} according to the new arguments given,
  then re-fit it.
  The actual re-fitting is performed by the model-fitting
  function \code{\link{ppm}}.

  If you are comparing several model fits to the same data,
  or fits of the same model to different data, it is
  strongly advisable to use \code{update.ppm}
  rather than trying to fit them by hand.
  This is because \code{update.ppm} re-fits the model
  in a way which is comparable to the original fit.

  The arguments \code{...} are matched to the formal arguments
  of \code{\link{ppm}} as follows.

  First, all the \emph{named} arguments in \code{...} are matched
  with the formal arguments of \code{\link{ppm}}.
  Use \code{name=NULL} to remove the argument \code{name} from the
  call.

  Second, any \emph{unnamed} arguments in \code{...} are
  matched with formal arguments of \code{\link{ppm}} if the matching
  is obvious from the class of the object. Thus \code{...} may contain
  \itemize{
    \item
    exactly one argument of class \code{"ppp"} or \code{"quad"},
    which will be interpreted as the named argument \code{Q};
    \item
    exactly one argument of class \code{"formula"}, which will be
    interpreted as the named argument \code{trend} (or as specifying
    a change to the trend formula);
    \item
    exactly one argument of class \code{"interact"}, which will be
    interpreted as the named argument \code{interaction};
    \item
    exactly one argument of class \code{"data.frame"}, which will be
    interpreted as the named argument \code{covariates}.
  }

  The \code{trend} argument can be a formula that specifies a
  \emph{change} to the current trend formula. For example, the
  formula \code{~ . + Z} specifies that the additional covariate
  \code{Z} will be added to the right hand side of the trend
  formula in the existing \code{object}.

  The argument \code{fixdummy=TRUE} ensures comparability of the
  objects before and after updating.
  When \code{fixdummy=FALSE}, calling \code{update.ppm}
  is exactly the same as calling \code{ppm} with the updated
  arguments. However, the original and updated models
  are not strictly comparable (for example, their pseudolikelihoods
  are not strictly comparable) unless they used the same set of dummy
  points for the quadrature scheme. Setting \code{fixdummy=TRUE}
  ensures that the re-fitting will be performed using the same set
  of dummy points. This is highly recommended.

  The value of \code{use.internal} determines where to find data
  to re-evaluate the model (data for the arguments mentioned in
  the original call to \code{ppm} that are not overwritten by
  arguments to \code{update.ppm}).
  
  If \code{use.internal=FALSE}, then arguments 
  to \code{ppm} are \emph{re-evaluated} in the frame where you
  call \code{update.ppm}. This is like the behaviour of the
  other methods for \code{\link{update}}. This means that if you have changed
  any of the objects referred to in the call, these changes will be
  taken into account. Also if the original call to \code{ppm} included
  any calls to random number generators, these calls will be recomputed,
  so that you will get a different outcome of the random numbers.

  If \code{use.internal=TRUE}, then arguments to \code{ppm} are extracted
  from internal data stored inside the current fitted
  model \code{object}. This is useful if you don't want to 
  re-evaluate anything. It is also necessary if 
  if \code{object} has been restored from a dump file
  using \code{\link{load}} or \code{\link{source}}. In such cases,
  we have lost the environment in which \code{object} was fitted,
  and data cannot be re-evaluated.

  By default, if \code{use.internal} is missing, \code{update.ppm} will
  re-evaluate the arguments if this is possible, and use internal data
  if not.
}
\value{
  Another fitted point process model (object of class \code{"ppm"}).
}
\examples{
  data(nztrees)
  data(cells)

  # fit the stationary Poisson process
  fit <- ppm(nztrees, ~ 1)

  # fit a nonstationary Poisson process
  fitP <- update(fit, trend=~x)
  fitP <- update(fit, ~x)

  # change the trend formula: add another term to the trend
  fitPxy <- update(fitP, ~ . + y)
  # change the trend formula: remove the x variable
  fitPy <- update(fitPxy, ~ . - x)

  # fit a stationary Strauss process
  fitS <- update(fit, interaction=Strauss(13))
  fitS <- update(fit, Strauss(13))

  # refit using a different edge correction
  fitS <- update(fitS, correction="isotropic")

  # re-fit the model to a subset
  # of the original point pattern
  nzw <- owin(c(0,148),c(0,95))
  nzsub <- nztrees[,nzw]
  fut <- update(fitS, Q=nzsub)
  fut <- update(fitS, nzsub)

  # WARNING: the point pattern argument is called 'Q'

  ranfit <- ppm(rpoispp(42), ~1, Poisson())
  ranfit
  # different random data!  
  update(ranfit)
  # the original data
  update(ranfit, use.internal=TRUE)  

}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{models}
