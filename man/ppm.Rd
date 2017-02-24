\name{ppm}
\alias{ppm}
\alias{ppm.formula}
\concept{point process model}
\concept{Poisson point process}
\concept{Gibbs point process}
\title{
  Fit Point Process Model to Data
}
\description{
  Fits a point process model to an observed point pattern.
}
\usage{
   ppm(Q, \dots)

   \method{ppm}{formula}(Q, interaction=NULL, \dots, data=NULL, subset)
}
\arguments{
  \item{Q}{
    A \code{formula} in the \R language describing the model
    to be fitted. 
  }
  \item{interaction}{
    An object of class \code{"interact"}
    describing the point process interaction
    structure, or a function that makes such an object,
    or \code{NULL} indicating that a Poisson process (stationary
    or nonstationary) should be fitted.
  }
  \item{\dots}{
    Arguments passed to \code{\link{ppm.ppp}}
    or \code{\link{ppm.quad}} to control the model-fitting process.
  }
  \item{data}{
    Optional. The values of spatial covariates (other than the Cartesian
    coordinates) required by the model.
    Either a data frame, or a list whose entries are images,
    functions, windows, tessellations or single numbers. See Details.
  }
  \item{subset}{
    Optional. An expression (typically involving the Cartesian coordinates
    \code{x} and \code{y} or the names of entries in \code{data})
    defining a subset of the spatial domain,
    to which the model-fitting should be restricted.
  }
}
\value{
  An object of class \code{"ppm"} describing a fitted point process
  model.
 
  See \code{\link{ppm.object}} for details of the format of this object
  and methods available for manipulating it.
}
\details{
  This function fits a point process model
  to an observed point pattern.
  The model may include
  spatial trend, interpoint interaction, and dependence on covariates.

  The model fitted by \code{ppm}
  is either a Poisson point process (in which different points
  do not interact with each other) or a Gibbs point process (in which
  different points typically inhibit each other).
  For clustered point process models, use \code{\link{kppm}}.

  The function \code{ppm} is generic, with methods for
  the classes \code{formula}, \code{ppp} and \code{quad}.
  This page describes the method for a \code{formula}.

  The first argument is a \code{formula} in the \R language
  describing the spatial trend model to be fitted. It has the general form
  \code{pattern ~ trend} where the left hand side \code{pattern} is usually
  the name of a spatial point pattern (object of class \code{"ppp"})
  to which the model should be fitted, or an expression which evaluates
  to a point pattern;
  and the right hand side \code{trend} is an expression specifying the
  spatial trend of the model.

  Systematic effects (spatial trend and/or dependence on 
  spatial covariates) are specified by the 
  \code{trend} expression on the right hand side of the formula.
  The trend may involve
  the Cartesian coordinates \code{x}, \code{y},
  the marks \code{marks},
  the names of entries in the argument \code{data} (if supplied),
  or the names of objects that exist in the \R session.
  The trend formula specifies the \bold{logarithm} of the
  intensity of a Poisson process, or in general, the logarithm of
  the first order potential of the Gibbs process.
  The formula should not use any names beginning with \code{.mpl}
  as these are reserved for internal use.
  If the formula is \code{pattern~1}, then
  the model to be fitted is stationary (or at least, its first order 
  potential is constant).
  
  The symbol \code{.} in the trend expression stands for
  all the covariates supplied in the argument \code{data}.
  For example the formula \code{pattern ~ .} indicates an additive
  model with a main effect for each covariate in \code{data}.
  
  Stochastic interactions between random points of the point process
  are defined by the argument \code{interaction}. This is an object of
  class \code{"interact"} which is initialised in a very similar way to the
  usage of family objects in \code{\link{glm}} and \code{gam}.
  The interaction models currently available are:
  \GibbsInteractionsList.
  See the examples below.
  Note that it is possible to combine several interactions
  using \code{\link{Hybrid}}.
 
  If \code{interaction} is missing or \code{NULL},
  then the model to be fitted
  has no interpoint interactions, that is, it is a Poisson process
  (stationary or nonstationary according to \code{trend}). In this case
  the methods of maximum pseudolikelihood and maximum logistic likelihood
  coincide with maximum likelihood. 
  
  The fitted point process model returned by this function can be printed 
  (by the print method \code{\link{print.ppm}})
  to inspect the fitted parameter values.
  If a nonparametric spatial trend was fitted, this can be extracted using
  the predict method \code{\link{predict.ppm}}.

  To fit a model involving spatial covariates
  other than the Cartesian coordinates \eqn{x} and \eqn{y},
  the values of the covariates should either be supplied in the
  argument \code{data}, or should be stored in objects that exist
  in the \R session.
  Note that it is not sufficient to have observed
  the covariate only at the points of the data point pattern; 
  the covariate must also have been observed at other 
  locations in the window.

  If it is given, the argument \code{data} is typically
  a list, with names corresponding to variables in the \code{trend} formula.
  Each entry in the list is either
  \describe{
    \item{a pixel image,}{
      giving the values of a spatial covariate at 
      a fine grid of locations. It should be an object of
      class \code{"im"}, see \code{\link{im.object}}.
    }
    \item{a function,}{
      which can be evaluated
      at any location \code{(x,y)} to obtain the value of the spatial
      covariate. It should be a \code{function(x, y)}
      or \code{function(x, y, ...)} in the \R language.
      The first two arguments of the function should be the
      Cartesian coordinates \eqn{x} and \eqn{y}. The function may have
      additional arguments; if the function does not have default
      values for these additional arguments, then the user must
      supply values for them, in \code{covfunargs}.
      See the Examples.
    }
    \item{a window,}{
      interpreted as a logical variable
      which is \code{TRUE} inside the window and \code{FALSE} outside
      it. This should be an object of class \code{"owin"}.
    }
    \item{a tessellation,}{
      interpreted as a factor covariate.
      For each spatial location, the factor value indicates
      which tile of the tessellation it belongs to.
      This should be an object of class \code{"tess"}.
    }
    \item{a single number,}{indicating a covariate that is
      constant in this dataset.
    }
  }
  The software will look up
  the values of each covariate at the required locations
  (quadrature points).

  Note that, for covariate functions, only the \emph{name} of the
  function appears in the trend formula. A covariate function is
  treated as if it were a single variable. The function arguments do not
  appear in the trend formula. See the Examples.

  If \code{data} is a list,
  the list entries should have names corresponding to
  (some of) the names of covariates in the model formula \code{trend}.
  The variable names \code{x}, \code{y} and \code{marks}
  are reserved for the Cartesian 
  coordinates and the mark values,
  and these should not be used for variables in \code{data}.

  Alternatively, \code{data} may be a data frame
  giving the values of the covariates at specified locations.
  Then \code{pattern} should be a quadrature scheme (object of class
  \code{"quad"}) giving the corresponding locations.
  See \code{\link{ppm.quad}} for details.
}
\section{Interaction parameters}{
  Apart from the Poisson model, every point process model fitted by
  \code{ppm} has parameters that determine the strength and
  range of \sQuote{interaction} or dependence between points.
  These parameters are of two types:
  \describe{
    \item{regular parameters:}{
      A parameter \eqn{\phi}{phi} is called \emph{regular}
      if the log likelihood is a linear function of \eqn{\theta}{theta} where 
      \eqn{\theta = \theta(\psi)}{theta = theta(psi)} is some transformation of 
      \eqn{\psi}{psi}. [Then \eqn{\theta}{theta} is called the canonical
      parameter.]
    }
    \item{irregular parameters}{
      Other parameters are called \emph{irregular}. 
    }
  }
  Typically, regular parameters determine the \sQuote{strength}
  of the interaction, while irregular parameters determine the
  \sQuote{range} of the interaction. For example, the Strauss process
  has a regular parameter \eqn{\gamma}{gamma} controlling the strength
  of interpoint inhibition, and an irregular parameter \eqn{r}
  determining the range of interaction.

  The \code{ppm} command is only designed to estimate regular
  parameters of the interaction.
  It requires the values of any irregular parameters of the interaction
  to be fixed. For example, to fit a Strauss process model to the \code{cells}
  dataset, you could type \code{ppm(cells ~ 1, Strauss(r=0.07))}.
  Note that the value of the irregular parameter \code{r} must be given.
  The result of this command will be a fitted model in which the
  regular parameter \eqn{\gamma}{gamma} has been estimated.

  To determine the irregular parameters, there are several
  practical techniques, but no general statistical theory available.
  Useful techniques include maximum profile pseudolikelihood, which
  is implemented in the command \code{\link{profilepl}},
  and Newton-Raphson maximisation, implemented in the
  experimental command \code{\link{ippm}}. 

  Some irregular parameters can be estimated directly from data:
  the hard-core radius in the model \code{\link{Hardcore}}
  and the matrix of hard-core radii in \code{\link{MultiHard}} can be
  estimated easily from data. In these cases, \code{ppm} allows the user
  to specify the interaction without giving
  the value of the irregular parameter. The user can give the
  hard core interaction as \code{interaction=Hardcore()}
  or even \code{interaction=Hardcore}, and 
  the hard core radius will then be estimated from the data.
}
\section{Technical Warnings and Error Messages}{
  See \code{\link{ppm.ppp}} for some technical warnings about the
  weaknesses of the algorithm, and explanation of some common error messages.
}
\references{
  Baddeley, A., Coeurjolly, J.-F., Rubak, E. and Waagepetersen, R. (2014)
  Logistic regression for spatial Gibbs point processes.
  \emph{Biometrika} \bold{101} (2) 377--392.

  Baddeley, A. and Turner, R. (2000)
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42} 283--322.
 
  Berman, M. and Turner, T.R. (1992)
  Approximating point process likelihoods with GLIM.
  \emph{Applied Statistics} \bold{41},  31--38.
 
  Besag, J. (1975)
  Statistical analysis of non-lattice data.
  \emph{The Statistician} \bold{24}, 179-195.
 
  Diggle, P.J., Fiksel, T., Grabarnik, P., Ogata, Y., Stoyan, D. and
  Tanemura, M. (1994)
  On parameter estimation for pairwise interaction processes.
  \emph{International Statistical Review} \bold{62}, 99-117.

  Huang, F. and Ogata, Y. (1999)
  Improvements of the maximum pseudo-likelihood estimators
  in various spatial statistical models.
  \emph{Journal of Computational and Graphical Statistics}
  \bold{8}, 510--530.
  
  Jensen, J.L. and Moeller, M. (1991)
  Pseudolikelihood for exponential family models of spatial point processes.
  \emph{Annals of Applied Probability} \bold{1}, 445--461.
 
  Jensen, J.L. and Kuensch, H.R. (1994)
  On asymptotic normality of pseudo likelihood
  estimates for pairwise interaction processes,
  \emph{Annals of the Institute of Statistical Mathematics}
  \bold{46}, 475--486.
}
\seealso{
  \code{\link{ppm.ppp}} and \code{\link{ppm.quad}} for
  more details on the fitting technique and edge correction.

  \code{\link{ppm.object}} for details of how to
  print, plot and manipulate a fitted model.

  \code{\link{ppp}} and \code{\link{quadscheme}}
  for constructing data.
  
  Interactions: \GibbsInteractionsList.

  See \code{\link{profilepl}} for advice on
  fitting nuisance parameters in the interaction,
  and \code{\link{ippm}} for irregular parameters in the trend.

  See \code{\link{valid.ppm}} and \code{\link{project.ppm}} for
  ensuring the fitted model is a valid point process.

  See \code{\link{kppm}} for fitting Cox point process models
  and cluster point process models, and \code{\link{dppm}} for fitting
  determinantal point process models.
}
\examples{
 # fit the stationary Poisson process
 # to point pattern 'nztrees'

 ppm(nztrees ~ 1)

 \dontrun{
 Q <- quadscheme(nztrees) 
 ppm(Q ~ 1) 
 # equivalent.
 }

fit1 <- ppm(nztrees ~ x)
 # fit the nonstationary Poisson process 
 # with intensity function lambda(x,y) = exp(a + bx)
 # where x,y are the Cartesian coordinates
 # and a,b are parameters to be estimated

fit1
coef(fit1)
coef(summary(fit1))

\dontrun{
 ppm(nztrees ~ polynom(x,2))
}
\testonly{
 ppm(nztrees ~ polynom(x,2), nd=16)
}

 # fit the nonstationary Poisson process 
 # with intensity function lambda(x,y) = exp(a + bx + cx^2)

 \dontrun{
 library(splines)
 ppm(nztrees ~ bs(x,df=3))
 }
 #       WARNING: do not use predict.ppm() on this result
 # Fits the nonstationary Poisson process 
 # with intensity function lambda(x,y) = exp(B(x))
 # where B is a B-spline with df = 3

\dontrun{
 ppm(nztrees ~ 1, Strauss(r=10), rbord=10)
}
\testonly{
 ppm(nztrees ~ 1, Strauss(r=10), rbord=10, nd=16)
}
 # Fit the stationary Strauss process with interaction range r=10
 # using the border method with margin rbord=10

\dontrun{
 ppm(nztrees ~ x, Strauss(13), correction="periodic")
}
\testonly{
 ppm(nztrees ~ x, Strauss(13), correction="periodic", nd=16)
}
 # Fit the nonstationary Strauss process with interaction range r=13
 # and exp(first order potential) =  activity = beta(x,y) = exp(a+bx)
 # using the periodic correction.


# Compare Maximum Pseudolikelihood, Huang-Ogata and Variational Bayes fits:
\dontrun{ppm(swedishpines ~ 1, Strauss(9))}

\dontrun{ppm(swedishpines ~ 1, Strauss(9), method="ho")}
\testonly{ppm(swedishpines ~ 1, Strauss(9), method="ho", nd=16, nsim=8)}

ppm(swedishpines ~ 1, Strauss(9), method="VBlogi")

 # COVARIATES
 #
 X <- rpoispp(42)
 weirdfunction <- function(x,y){ 10 * x^2 + 5 * sin(10 * y) }
 #
 # (a) covariate values as function
 ppm(X ~ y + weirdfunction)
 #
 # (b) covariate values in pixel image
 Zimage <- as.im(weirdfunction, unit.square())
 ppm(X ~ y + Z, covariates=list(Z=Zimage))
 #
 # (c) covariate values in data frame
 Q <- quadscheme(X)
 xQ <- x.quad(Q)
 yQ <- y.quad(Q)
 Zvalues <- weirdfunction(xQ,yQ)
 ppm(Q ~  y + Z, data=data.frame(Z=Zvalues))
 # Note Q not X

 # COVARIATE FUNCTION WITH EXTRA ARGUMENTS
 #
f <- function(x,y,a){ y - a }
ppm(X ~ x + f, covfunargs=list(a=1/2))

 # COVARIATE: inside/outside window
 b <- owin(c(0.1, 0.6), c(0.1, 0.9))
 ppm(X ~ b)

 ## MULTITYPE POINT PROCESSES ### 
 # fit stationary marked Poisson process
 # with different intensity for each species
\dontrun{ppm(lansing ~  marks, Poisson())}
\testonly{
  ama <- amacrine[square(0.7)]
  a <- ppm(ama ~  marks, Poisson(), nd=16)
}

 # fit nonstationary marked Poisson process
 # with different log-cubic trend for each species
\dontrun{ppm(lansing ~  marks * polynom(x,y,3), Poisson())}
\testonly{b <- ppm(ama ~  marks * polynom(x,y,2), Poisson(), nd=16)}

}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{models}
