\name{addvar}
\alias{addvar}
\title{
  Added Variable Plot for Point Process Model
}
\description{
  Computes the coordinates for an Added Variable Plot
  for a fitted point process model.
}
\usage{
addvar(model, covariate, ...,
                   subregion=NULL,
                   bw="nrd0", adjust=1,
                   from=NULL, to=NULL, n=512,
                   bw.input = c("points", "quad"),
                   bw.restrict = FALSE,
                   covname, crosscheck=FALSE)
}
\arguments{
  \item{model}{
    Fitted point process model (object of class \code{"ppm"}).
  }
  \item{covariate}{
    The covariate to be added to the model. Either a
    pixel image, a \code{function(x,y)}, or a character string
    giving the name of a covariate that was supplied when
    the model was fitted.
  }
  \item{subregion}{
    Optional.  A window (object of class \code{"owin"})
    specifying a subset of the spatial domain of the data.
    The calculation will be confined to the data in this subregion.
  }
  \item{bw}{
    Smoothing bandwidth or bandwidth rule
    (passed to \code{\link[stats]{density.default}}).
  }
  \item{adjust}{
    Smoothing bandwidth adjustment factor
    (passed to \code{\link[stats]{density.default}}).
  }
  \item{n, from, to}{
    Arguments passed to \code{\link[stats]{density.default}} to
    control the number and range of values at which the function
    will be estimated.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link[stats]{density.default}}.
  }
  \item{bw.input}{
    Character string specifying the input data used for automatic
    bandwidth selection. 
  }
  \item{bw.restrict}{
    Logical value, specifying whether bandwidth selection is performed using
    data from the entire spatial domain or from the \code{subregion}.
  }
  \item{covname}{
    Optional. Character string to use as the name of the covariate.
  }
  \item{crosscheck}{
    For developers only.
    Logical value indicating whether to perform
    cross-checks on the validity of the calculation.
  }
}
\details{
  This command generates the plot coordinates for an Added Variable Plot
  for a spatial point process model.
  
  Added Variable Plots (Cox, 1958, sec 4.5; Wang, 1985)
  are commonly used in linear models and generalized linear
  models, to decide whether a model with response \eqn{y} and predictors \eqn{x}
  would be improved by including another predictor \eqn{z}.
  
  In a (generalised) linear model 
  with response \eqn{y} and predictors \eqn{x},
  the Added Variable Plot for a new covariate \eqn{z} 
  is a plot of the smoothed Pearson residuals from the original model
  against the scaled residuals from a weighted linear
  regression of \eqn{z} on \eqn{x}.
  If this plot has nonzero slope, then the new covariate \eqn{z} is
  needed. For general advice see Cook and Weisberg(1999); Harrell (2001).
  
  Essentially the same technique can be used for a spatial point process
  model (Baddeley et al, 2012).

  The argument \code{model} should be a fitted spatial point process
  model (object of class \code{"ppm"}). 

  The argument \code{covariate}
  identifies the covariate that is to be considered for addition to
  the model. It should be either a pixel image (object of class
  \code{"im"}) or a \code{function(x,y)} giving the values of the
  covariate at any spatial location. Alternatively \code{covariate}
  may be a character string, giving the name of a covariate that was
  supplied (in the \code{covariates} argument to \code{\link{ppm}})
  when the model was fitted, but was not used in the model.

  The result of \code{addvar(model, covariate)} is an object belonging
  to the classes \code{"addvar"} and \code{"fv"}. Plot this object to
  generate the added variable plot. 
  
  Note that the plot method shows the pointwise significance bands
  for a test of the \emph{null} model, i.e. the null hypothesis
  that the new covariate has no effect.

  The smoothing bandwidth is controlled by the arguments
  \code{bw}, \code{adjust}, \code{bw.input} and \code{bw.restrict}.
  If \code{bw} is a numeric value, then
  the bandwidth is taken to be \code{adjust * bw}.
  If \code{bw} is a string representing a bandwidth selection rule
  (recognised by \code{\link[stats]{density.default}})
  then the bandwidth is selected by this rule.

  The data used for automatic bandwidth selection are
  specified by \code{bw.input} and \code{bw.restrict}.
  If \code{bw.input="points"}  (the default) then bandwidth selection is
  based on the covariate values at the points of the original point
  pattern dataset to which the model was fitted.
  If \code{bw.input="quad"} then bandwidth selection is
  based on the covariate values at every quadrature point used to
  fit the model.
  If \code{bw.restrict=TRUE} then the bandwidth selection is performed
  using only data from inside the \code{subregion}.
}
\section{Slow computation}{
  In a large dataset, computation can be very slow if the default
  settings are used, because the smoothing bandwidth is selected
  automatically. To avoid this, specify a numerical value
  for the bandwidth \code{bw}. One strategy is to use a coarser
  subset of the data to select \code{bw} automatically.
  The selected bandwidth can be read off the print output for
  \code{addvar}.  
}
\value{
  An object of class \code{"addvar"} containing the coordinates
  for the added variable plot. There is a \code{plot} method.
}
\section{Internal data}{
  The return value has an attribute \code{"spatial"} which contains
  the internal data: the computed values of the residuals,
  and of all relevant covariates,
  at each quadrature point of the model. It is an object of class
  \code{"ppp"} with a data frame of marks.
}
\references{
  Baddeley, A., Chang, Y.-M., Song, Y. and Turner, R. (2013)
  Residual diagnostics for covariate effects in
  spatial point process models.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{22}, 886--905.
  
  Cook, R.D. and Weisberg, S. (1999)
  \emph{Applied regression, including computing and graphics}.
  New York: Wiley.
  
  Cox, D.R. (1958) \emph{Planning of Experiments}. New York: Wiley.

  Harrell, F. (2001) \emph{Regression Modeling Strategies}. New York: Springer.

  Wang, P. (1985) Adding a variable in generalized linear models.
  \emph{Technometrics} \bold{27}, 273--276.
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  Rolf Turner
  \email{r.turner@auckland.ac.nz},
  Ya-Mei Chang and Yong Song.
}
\seealso{
  \code{\link{parres}},
  \code{\link{rhohat}},
  \code{\link{rho2hat}}.
}
\examples{
  X <-  rpoispp(function(x,y){exp(3+3*x)})
  model <- ppm(X, ~y)
  adv <- addvar(model, "x")
  plot(adv)
  adv <- addvar(model, "x", subregion=square(0.5))
}
\keyword{spatial}
\keyword{models}
