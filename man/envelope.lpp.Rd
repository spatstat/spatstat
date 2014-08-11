\name{envelope.lpp}
\alias{envelope.lpp}
\alias{envelope.lppm}
\title{
  Envelope for Point Patterns on Linear Network
}
\description{
  Enables envelopes to be computed for point patterns on a linear network.
}
\usage{
  \method{envelope}{lpp}(Y, fun=linearK, nsim=99, nrank=1, \dots, 
  simulate=NULL, verbose=TRUE, 
  transform=NULL,global=FALSE,ginterval=NULL,
  savefuns=FALSE, savepatterns=FALSE,
  nsim2=nsim, VARIANCE=FALSE, nSD=2, Yname=NULL,
  do.pwrong=FALSE, envir.simul=NULL)


  \method{envelope}{lppm}(Y, fun=linearK, nsim=99, nrank=1, \dots, 
  simulate=NULL, verbose=TRUE, 
  transform=NULL,global=FALSE,ginterval=NULL,
  savefuns=FALSE, savepatterns=FALSE,
  nsim2=nsim, VARIANCE=FALSE, nSD=2, Yname=NULL,
  do.pwrong=FALSE, envir.simul=NULL)
}
\arguments{
  \item{Y}{
    A point pattern on a linear network
    (object of class \code{"lpp"})
    or a fitted point process model on a linear network
    (object of class \code{"lppm"}).
  }
  \item{fun}{
    Function that is to be computed for each simulated pattern.
  }
  \item{nsim}{
    Number of simulations to perform.
  }
  \item{nrank}{
    Integer. Rank of the envelope value amongst the \code{nsim} simulated
    values. A rank of 1 means that the minimum and maximum
    simulated values will be used.
  }
   \item{\dots}{
    Extra arguments passed to \code{fun}.
  }
  \item{simulate}{
    Optional. Specifies how to generate the simulated point patterns.
    If \code{simulate} is an expression in the R language, then this
    expression will be evaluated \code{nsim} times,
    to obtain \code{nsim} point patterns which are taken as the
    simulated patterns from which the envelopes are computed.
    If \code{simulate} is a list of point patterns, then the entries
    in this list will be treated as the simulated patterns from which
    the envelopes are computed.
    Alternatively \code{simulate} may be an object produced by the
    \code{envelope} command: see Details.
  }
  \item{verbose}{
    Logical flag indicating whether to print progress reports
    during the simulations.
  }
  \item{transform}{
    Optional. A transformation to be applied to the
    function values, before the envelopes are computed.
    An expression object (see Details).
  }
  \item{global}{
    Logical flag indicating whether envelopes should be pointwise
    (\code{global=FALSE}) or simultaneous (\code{global=TRUE}).
  }
  \item{ginterval}{
    Optional.
    A vector of length 2 specifying
    the interval of \eqn{r} values for the simultaneous critical
    envelopes. Only relevant if \code{global=TRUE}.
  }
  \item{savefuns}{
    Logical flag indicating whether to save all the simulated
    function values.
  }
  \item{savepatterns}{
    Logical flag indicating whether to save all the simulated
    point patterns.
  }
  \item{nsim2}{
    Number of extra simulated point patterns to be generated
    if it is necessary to use simulation to estimate the theoretical
    mean of the summary function. Only relevant when \code{global=TRUE}
    and the simulations are not based on CSR.
  }
  \item{VARIANCE}{
    Logical. If \code{TRUE}, critical envelopes will be calculated
    as sample mean plus or minus \code{nSD} times sample standard
    deviation.
  }
  \item{nSD}{
    Number of estimated standard deviations used to determine
    the critical envelopes, if \code{VARIANCE=TRUE}.
  }
  \item{Yname}{
    Character string that should be used as the name of the 
    data point pattern \code{Y} when printing or plotting the results.
  }
  \item{do.pwrong}{
    Logical. If \code{TRUE}, the algorithm will also estimate
    the true significance level of the \dQuote{wrong} test (the test that
    declares the summary function for the data to be significant
    if it lies outside the \emph{pointwise} critical boundary at any
    point). This estimate is printed when the result is printed.
  }
  \item{envir.simul}{
    Environment in which to evaluate the expression \code{simulate},
    if not the current environment.
  }
}
\details{
  This is a method for the generic
  function \code{\link[spatstat]{envelope}} 
  applicable to point patterns on a linear network.
  
  The argument \code{Y} can be either a point pattern on a linear
  network, or a fitted point process model on a linear network. 
  The function \code{fun} will be evaluated for the data
  and also for \code{nsim} simulated point
  patterns on the same linear network.
  The upper and lower
  envelopes of these evaluated functions will be computed
  as described in \code{\link[spatstat]{envelope}}.
  
  The type of simulation is determined as follows.
  \itemize{
    \item
    if \code{Y} is a point pattern (object of class \code{"lpp"})
    and \code{simulate} is missing or \code{NULL},
    then random point patterns will be generated according to
    a Poisson point process on the linear network on which \code{Y}
    is defined, with intensity estimated from \code{Y}.
    \item
    if \code{Y} is a fitted point process model (object of class
    \code{"lppm"}) and \code{simulate} is missing or \code{NULL},
    then random point patterns will be generated by simulating
    from the fitted model.
    \item 
    If \code{simulate} is present, it should be an expression that
    can be evaluated to yield random point patterns on the same
    linear network as \code{Y}.
  }
  
  The function \code{fun} should accept as its first argument
  a point pattern on a linear network (object of class \code{"lpp"})
  and should have another argument called \code{r} or a \code{\dots}
  argument. 
}
\value{
  Function value table (object of class \code{"fv"})
  with additional information,
  as described in \code{\link[spatstat]{envelope}}.
}
\author{
   Ang Qi Wei \email{aqw07398@hotmail.com} and
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
   \code{\link[spatstat]{envelope}},
   \code{\link{linearK}}
}
\references{
  Ang, Q.W. (2010)
  \emph{Statistical methodology for events on a network}.
  Master's thesis, School of Mathematics and Statistics, University of
  Western Australia.
  
  Ang, Q.W., Baddeley, A. and Nair, G. (2012)
  Geometrically corrected second-order analysis of 
  events on a linear network, with applications to
  ecology and criminology.
  To appear in \emph{Scandinavian Journal of Statistics}.
  
  Okabe, A. and Yamada, I. (2001) The K-function method on a network and
  its computational implementation. \emph{Geographical Analysis}
  \bold{33}, 271-290.
}
\examples{
   if(interactive()) {
     ns <- 39
     np <- 40
   } else { ns <- np <- 3 }
   X <- runiflpp(np, simplenet)
   # uniform Poisson
   envelope(X, nsim=ns)
   # nonuniform Poisson
   fit <- lppm(X, ~x)
   envelope(fit, nsim=ns)

   #multitype
   marks(X) <- sample(letters[1:2], np, replace=TRUE)
   envelope(X, nsim=ns)
}
\keyword{spatial}
