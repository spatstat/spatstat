\name{profilepl}
\alias{profilepl}
\title{Profile Maximum Pseudolikelihood}
\description{
  Fits point process models by profile maximum pseudolikelihood
}
\usage{
profilepl(s, f, \dots, rbord = NULL, verbose = TRUE)
}
\arguments{
  \item{s}{
    Data frame containing values of the irregular parameters
    over which the profile pseudolikelihood will be computed.
  }
  \item{f}{
    Function (such as \code{\link{Strauss}})
    that generates an interpoint interaction object, given
    values of the irregular parameters.
  }
  \item{\dots}{
    Data passed to \code{\link{ppm}} to fit the model.
  }
  \item{rbord}{
    Radius for border correction (same for all models).
    If omitted, this will be computed from the interactions.
  }
  \item{verbose}{
    Logical flag indicating whether to print progress reports.
  }
}
\details{
  The model-fitting function \code{\link{ppm}} fits point process
  models to point pattern data. However, 
  only the \sQuote{regular} parameters of the model can be fitted by
  \code{\link{ppm}}. The model may also depend on \sQuote{irregular}
  parameters that must be fixed in any call to \code{\link{ppm}}.

  This function \code{profilepl} is a wrapper which finds the values of the
  irregular parameters that give the best fit. It uses the method of
  maximum profile pseudolikelihood. 
  
  The argument \code{s} must be a data frame whose columns contain
  values of the irregular parameters over which the maximisation is
  to be performed.

  An irregular parameter may affect either the interpoint interaction
  or the spatial trend. 
  
  \describe{
    \item{interaction parameters:}{
      in a call to \code{\link{ppm}}, the argument \code{interaction}
      determines the interaction between points. It is usually
      a call to a function such as \code{\link{Strauss}}. The
      arguments of this call are irregular parameters.
      For example, the interaction radius parameter \eqn{r} of the Strauss
      process, determined by the argument \code{r}
      to the function \code{\link{Strauss}}, is an irregular parameter.
    }
    \item{trend parameters:}{
      in a call to \code{\link{ppm}}, the spatial trend may depend on
      covariates, which are supplied by the argument \code{covariates}.
      These covariates may be functions written by the user,
      of the form \code{function(x,y,...)}, and the extra arguments
      \code{\dots} are irregular parameters.
    }
  }
      
  The argument \code{f} determines the interaction
  for each model to be fitted. It would typically be one of the functions
  \code{\link{Poisson}},
  \code{\link{AreaInter}},
  \code{\link{BadGey}},
  \code{\link{DiggleGatesStibbard}},
  \code{\link{DiggleGratton}},
  \code{\link{Fiksel}},
  \code{\link{Geyer}},
  \code{\link{Hardcore}},
  \code{\link{LennardJones}},
  \code{\link{OrdThresh}}, 
  \code{\link{Softcore}},
  \code{\link{Strauss}} or
  \code{\link{StraussHard}}.
  Alternatively it could be a function written by the user.

  Columns of \code{s} which match the names of arguments of \code{f}
  will be interpreted as interaction parameters. Other columns will be
  interpreted as trend parameters.

  To apply the method of profile maximum pseudolikelihood,
  each row of \code{s} will be taken in turn. Interaction parameters in this
  row will be passed to \code{f}, resulting in an interaction object. 
  Then \code{\link{ppm}} will be applied to the data \code{...}
  using this interaction. Any trend parameters will be passed to
  \code{\link{ppm}} through the argument \code{covfunargs}.
  This results in a fitted point process model.
  The value of the log pseudolikelihood from this model is stored.
  After all rows of \code{s} have been processed in this way, the
  row giving the maximum value of log pseudolikelihood will be found.

  The object returned by \code{profilepl} contains the profile
  pseudolikelihood function, the best fitting model, and other data.
  It can be plotted (yielding a
  plot of the log pseudolikelihood values against the irregular
  parameters) or printed (yielding information about the best fitting
  values of the irregular parameters). 
 
  In general, \code{f} may be any function that will return
  an interaction object (object of class \code{"interact"})
  that can be used in a call to \code{\link{ppm}}. Each argument of
  \code{f} must be a single value. 
}
\value{
  An object of class \code{"profilepl"}. There are methods
  for \code{\link{plot}} and \code{\link{print}} for this class.

  The components of the object include
  \item{fit}{Best-fitting model}
  \item{param}{The data frame \code{s}}
  \item{iopt}{Row index of the best-fitting parameters in \code{s}}

  To extract the best fitting model you can also use \code{\link{as.ppm}}.
}
\examples{
    # one irregular parameter
    s <- data.frame(r=seq(0.05,0.15, by=0.01))
    \testonly{
      s <- data.frame(r=c(0.05,0.1,0.15))
    }
    ps <- profilepl(s, Strauss, cells)
    ps
    if(interactive()) plot(ps)

    # two irregular parameters
    s <- expand.grid(r=seq(0.05,0.15, by=0.01),sat=1:3)
    \testonly{
      s <- expand.grid(r=c(0.07,0.12),sat=1:2)
    }
    pg <- profilepl(s, Geyer, cells)
    pg
    if(interactive()) {
      plot(pg)
      as.ppm(fit)
    }

    # multitype pattern with a common interaction radius
    \dontrun{
     s <- data.frame(R=seq(0.03,0.05,by=0.01))
     MS <- function(R) { MultiStrauss(radii=diag(c(R,R))) }
     pm <- profilepl(s, MS, amacrine, ~marks)
    }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}

