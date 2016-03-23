\name{profilepl}
\alias{profilepl}
\title{Fit Models by Profile Maximum Pseudolikelihood or AIC}
\description{
  Fits point process models by maximising the profile likelihood,
  profile pseudolikelihood, profile composite likelihood or AIC.
}
\usage{
profilepl(s, f, \dots, aic=FALSE, rbord=NULL, verbose = TRUE)
}
\arguments{
  \item{s}{
    Data frame containing values of the irregular parameters
    over which the criterion will be computed.
  }
  \item{f}{
    Function (such as \code{\link{Strauss}})
    that generates an interpoint interaction object, given
    values of the irregular parameters.
  }
  \item{\dots}{
    Data passed to \code{\link{ppm}} to fit the model.
  }
  \item{aic}{
    Logical value indicating whether to find the parameter values
    which minimise the AIC (\code{aic=TRUE}) or maximise the
    profile likelihood (\code{aic=FALSE}, the default).
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
  irregular parameters that give the best fit.
  If \code{aic=FALSE} (the default),
  the best fit is the model which maximises the
  likelihood (if the models are Poisson processes) or maximises
  the pseudolikelihood or logistic likelihood. 
  If \code{aic=TRUE} then the best fit is the model which
  minimises the Akaike Information Criterion \code{\link{AIC.ppm}}.
  
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

  The data frame \code{s} must provide values for each argument of
  \code{f}, except for the optional arguments, which are those arguments of
  \code{f} that have the default value \code{NA}.  

  To find the best fit,
  each row of \code{s} will be taken in turn. Interaction parameters in this
  row will be passed to \code{f}, resulting in an interaction object. 
  Then \code{\link{ppm}} will be applied to the data \code{...}
  using this interaction. Any trend parameters will be passed to
  \code{\link{ppm}} through the argument \code{covfunargs}.
  This results in a fitted point process model.
  The value of the log pseudolikelihood or AIC from this model is stored.
  After all rows of \code{s} have been processed in this way, the
  row giving the maximum value of log pseudolikelihood will be found.

  The object returned by \code{profilepl} contains the profile
  pseudolikelihood (or profile AIC) function,
  the best fitting model, and other data.
  It can be plotted (yielding a
  plot of the log pseudolikelihood or AIC values against the irregular
  parameters) or printed (yielding information about the best fitting
  values of the irregular parameters). 
 
  In general, \code{f} may be any function that will return
  an interaction object (object of class \code{"interact"})
  that can be used in a call to \code{\link{ppm}}. Each argument of
  \code{f} must be a single value. 
}
\value{
  An object of class \code{"profilepl"}. There are methods
  for \code{\link[graphics]{plot}},
  \code{\link[base]{print}},
  \code{\link[base]{summary}},
  \code{\link[stats]{simulate}},
  \code{\link{as.ppm}}
  and 
  \code{\link{parameters}} for objects of this class.

  The components of the object include
  \item{fit}{Best-fitting model}
  \item{param}{The data frame \code{s}}
  \item{iopt}{Row index of the best-fitting parameters in \code{s}}

  To extract the best fitting model you can also use \code{\link{as.ppm}}.
}
\examples{
    # one irregular parameter
    rr <- data.frame(r=seq(0.05,0.15, by=0.01))
    \testonly{
      rr <- data.frame(r=c(0.05,0.1,0.15))
    }
    ps <- profilepl(rr, Strauss, cells)
    ps
    if(interactive()) plot(ps)

    # two irregular parameters
    rs <- expand.grid(r=seq(0.05,0.15, by=0.01),sat=1:3)
    \testonly{
      rs <- expand.grid(r=c(0.07,0.12),sat=1:2)
    }
    pg <- profilepl(rs, Geyer, cells)
    pg
    if(interactive()) {
      plot(pg)
      as.ppm(pg)
    }

    # multitype pattern with a common interaction radius
    \dontrun{
     RR <- data.frame(R=seq(0.03,0.05,by=0.01))
     MS <- function(R) { MultiStrauss(radii=diag(c(R,R))) }
     pm <- profilepl(RR, MS, amacrine ~marks)
    }
    ## more information 
    summary(pg)
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au},
 Rolf Turner \email{r.turner@auckland.ac.nz}
 and Ege Rubak \email{rubak@math.aau.dk}
}
\keyword{spatial}
\keyword{models}

