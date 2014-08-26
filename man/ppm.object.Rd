\name{ppm.object}
\alias{ppm.object}  %DoNotExport
\alias{methods.ppm} %DoNotExport
\title{Class of Fitted Point Process Models}
\description{
  A class \code{ppm} to represent a fitted stochastic model
  for a point process. The output of \code{\link{ppm}}.
}
\details{
  An object of class \code{ppm} represents a stochastic point process
  model that has been fitted to a point pattern dataset.
  Typically it is the output of the model fitter,
  \code{\link{ppm}}.

  The class \code{ppm} has methods for the following
  standard generic functions:

  \tabular{lll}{
    generic \tab method \tab description \cr
    \code{print} \tab \code{\link{print.ppm}}
    \tab print details \cr
    \code{plot} \tab \code{\link{plot.ppm}}
    \tab plot fitted model \cr
    \code{predict} \tab \code{\link{predict.ppm}}
    \tab fitted intensity and conditional intensity \cr 
    \code{fitted} \tab \code{\link{fitted.ppm}}
    \tab fitted intensity \cr
    \code{coef} \tab \code{\link{coef.ppm}}
    \tab fitted coefficients of model \cr
    \code{anova} \tab \code{\link{anova.ppm}}
    \tab Analysis of Deviance \cr
    \code{formula} \tab \code{\link{formula.ppm}}
    \tab Extract model formula \cr
    \code{terms} \tab \code{\link{terms.ppm}}
    \tab Terms in the model formula \cr
    \code{labels} \tab \code{labels.ppm}
    \tab Names of estimable terms in the model formula \cr
    \code{residuals} \tab \code{\link{residuals.ppm}}
    \tab Point process residuals \cr
    \code{simulate} \tab \code{\link{simulate.ppm}}
    \tab Simulate the fitted model \cr
    \code{update} \tab \code{\link{update.ppm}}
    \tab Change or refit the model \cr
    \code{vcov} \tab \code{\link{vcov.ppm}}
    \tab Variance/covariance matrix of parameter estimates \cr
    \code{model.frame} \tab \code{\link{model.frame.ppm}}
    \tab Model frame \cr
    \code{model.matrix} \tab \code{\link{model.matrix.ppm}}
    \tab Design matrix \cr
    \code{logLik} \tab \code{\link{logLik.ppm}}
    \tab log \emph{pseudo} likelihood \cr
    \code{extractAIC} \tab \code{\link{extractAIC.ppm}}
    \tab pseudolikelihood counterpart of AIC \cr
    \code{nobs} \tab \code{\link{nobs.ppm}}
    \tab number of observations 
  }

  Objects of class \code{ppm} can also be handled by the
  following standard functions, without requiring a special method:

  \tabular{ll}{
    name \tab description \cr
    \code{\link{confint}} \tab Confidence intervals for parameters \cr
    \code{\link{step}} \tab Stepwise model selection \cr
    \code{\link{drop1}} \tab One-step model improvement \cr
    \code{\link{add1}} \tab  One-step model improvement
  }

  The class \code{ppm} also has methods for the following
  generic functions defined in the \pkg{spatstat} package:

  \tabular{lll}{
    generic \tab method \tab description \cr
    \code{\link{as.interact}} \tab \code{\link{as.interact.ppm}}
    \tab Interpoint interaction structure \cr
    \code{\link{as.owin}} \tab \code{\link{as.owin.ppm}}
    \tab Observation window of data \cr
    \code{\link{berman.test}} \tab \code{\link{berman.test.ppm}}
    \tab Berman's test \cr
    \code{\link{envelope}} \tab \code{\link{envelope.ppm}}
    \tab Simulation envelopes \cr
    \code{\link{fitin}} \tab \code{\link{fitin.ppm}}
    \tab Fitted interaction \cr
    \code{\link{is.marked}} \tab \code{\link{is.marked.ppm}}
    \tab Determine whether the model is marked \cr
    \code{\link{is.multitype}} \tab \code{\link{is.multitype.ppm}}
    \tab Determine whether the model is multitype \cr
    \code{\link{is.poisson}} \tab \code{\link{is.poisson.ppm}}
    \tab Determine whether the model is Poisson \cr
    \code{\link{is.stationary}} \tab \code{\link{is.stationary.ppm}}
    \tab Determine whether the model is stationary \cr
    \code{\link{cdf.test}} \tab \code{\link{cdf.test.ppm}}
    \tab Spatial distribution test \cr
    \code{\link{quadrat.test}} \tab \code{\link{quadrat.test.ppm}}
    \tab Quadrat counting test \cr
    \code{\link{reach}} \tab \code{\link{reach.ppm}}
    \tab Interaction range of model \cr
    \code{\link{rmhmodel}} \tab \code{\link{rmhmodel.ppm}}
    \tab Model in a form that can be simulated \cr
    \code{\link{rmh}} \tab \code{\link{rmh.ppm}}
    \tab Perform simulation \cr
    \code{\link{unitname}} \tab \code{\link{unitname.ppm}}
    \tab Name of unit of length
    }
    
  Information about the data (to which the model was fitted)
  can be extracted using \code{\link{data.ppm}}, \code{\link{dummy.ppm}}
  and \code{\link{quad.ppm}}.
}
\section{Internal format}{
  If you really need to get at the internals,
  a \code{ppm} object contains at least the following entries:
  \tabular{ll}{
    \code{coef} \tab  the fitted regular parameters (as returned by
                     \code{glm}) \cr
    \code{trend} \tab  the trend formula or \code{NULL} \cr
    \code{interaction} \tab the point process interaction family 
                        (an object of class \code{"interact"})
			or \code{NULL} \cr
    \code{Q} \tab      the quadrature scheme used \cr
    \code{maxlogpl} \tab the maximised value of log pseudolikelihood \cr
    \code{correction} \tab  name of edge correction method used \cr
  }
  See \code{\link{ppm}} for explanation of these concepts.
  The irregular parameters (e.g. the interaction radius of the
  Strauss process) are encoded in the \code{interaction} entry.
  However see the Warnings.
}

\seealso{
  \code{\link{ppm}},
  \code{\link{coef.ppm}},
  \code{\link{fitted.ppm}},
  \code{\link{print.ppm}},
  \code{\link{predict.ppm}},
  \code{\link{plot.ppm}}.
}
\section{Warnings}{
  The internal representation of \code{ppm} objects
  may change slightly between releases of the \pkg{spatstat} package.
}
\examples{
  data(cells)
  fit <- ppm(cells, ~ x, Strauss(0.1), correction="periodic")
  fit
  coef(fit)
  \dontrun{
  pred <- predict(fit)
  }
  pred <- predict(fit, ngrid=20, type="trend")
  \dontrun{
  plot(fit)
  }  
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{attribute}
 
 
