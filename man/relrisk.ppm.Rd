\name{relrisk.ppm}
\alias{relrisk.ppm}
\title{
  Parametric Estimate of Spatially-Varying Relative Risk
}
\description{
  Given a point process model fitted to a multitype point pattern,
  this function computes the fitted spatially-varying probability
  of each type of point, or the ratios of such probabilities,
  according to the fitted model. Optionally the standard errors of the estimates
  are also computed.
}
\usage{
\method{relrisk}{ppm}(X, \dots,
                     at = c("pixels", "points"),
                     relative = FALSE, se = FALSE,
                     casecontrol = TRUE, control = 1, case,
                     ngrid = NULL, window = NULL)
}
\arguments{
  \item{X}{
    A fitted point process model (object of class \code{"ppm"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{at}{
    String specifying whether to compute the probability values
    at a grid of pixel locations (\code{at="pixels"}) or
    only at the points of \code{X} (\code{at="points"}).
  }
  \item{relative}{
    Logical.
    If \code{FALSE} (the default) the algorithm
    computes the probabilities of each type of point.
    If \code{TRUE}, it computes the    
    \emph{relative risk}, the ratio of probabilities
    of each type relative to the probability of a control.
  }
  \item{se}{
    Logical value indicating whether to compute standard errors
    as well.
  }
  \item{casecontrol}{
    Logical. Whether to treat a bivariate point pattern
    as consisting of cases and controls, and return only the
    probability or relative risk of a case.
    Ignored if there are more than 2 types of points.
    See Details.
  }
  \item{control}{
    Integer, or character string, identifying which mark value
    corresponds to a control. 
  }
  \item{case}{
    Integer, or character string, identifying which mark value
    corresponds to a case (rather than a control)
    in a bivariate point pattern.
    This is an alternative to the argument \code{control}
    in a bivariate point pattern. 
    Ignored if there are more than 2 types of points.
  }
  \item{ngrid}{
    Optional. Dimensions of a rectangular grid of locations
    inside \code{window} where the predictions should be computed.
    An integer, or an integer vector of length 2,
    specifying the number of grid points in the \eqn{y} and \eqn{x}
    directions.
    (Applies only when \code{at="pixels"}.)
  }
  \item{window}{
    Optional. A window (object of class \code{"owin"})
    \emph{delimiting} the locations where predictions
    should be computed. Defaults to the window of the
    original data used to fit the model \code{object}.
    (Applies only when \code{at="pixels"}.)
  }
}
\details{
  The command \code{\link{relrisk}} is generic and can be used to
  estimate relative risk in different ways.
  
  This function \code{relrisk.ppm} is the method for fitted point
  process models (class \code{"ppm"}). It computes \emph{parametric}
  estimates of relative risk, using the fitted model.

  If \code{X}  is a bivariate point pattern
  (a multitype point pattern consisting of two types of points)
  then by default,
  the points of the first type (the first level of \code{marks(X)})
  are treated as controls or non-events, and points of the second type
  are treated as cases or events. Then by default this command computes
  the spatially-varying \emph{probability} of a case,
  i.e. the probability \eqn{p(u)}
  that a point at spatial location \eqn{u}
  will be a case. If \code{relative=TRUE}, it computes the
  spatially-varying \emph{relative risk} of a case relative to a
  control, \eqn{r(u) = p(u)/(1- p(u))}.

  If \code{X} is a multitype point pattern with \eqn{m > 2} types,
  or if \code{X} is a bivariate point pattern
  and \code{casecontrol=FALSE},
  then by default this command computes, for each type \eqn{j},
  a nonparametric estimate of
  the spatially-varying \emph{probability} of an event of type \eqn{j}.
  This is the probability \eqn{p_j(u)}{p[j](u)}
  that a point at spatial location \eqn{u}
  will belong to type \eqn{j}.
  If \code{relative=TRUE}, the command computes the
  \emph{relative risk} of an event of type \eqn{j}
  relative to a control,
  \eqn{r_j(u) = p_j(u)/p_k(u)}{r[j](u) = p[j](u)/p[k](u)},
  where events of type \eqn{k} are treated as controls.
  The argument \code{control} determines which type \eqn{k}
  is treated as a control.

  If \code{at = "pixels"} the calculation is performed for
  every spatial location \eqn{u} on a fine pixel grid, and the result
  is a pixel image representing the function \eqn{p(u)}
  or a list of pixel images representing the functions 
  \eqn{p_j(u)}{p[j](u)} or \eqn{r_j(u)}{r[j](u)}
  for \eqn{j = 1,\ldots,m}{j = 1,...,m}.
  An infinite value of relative risk (arising because the
  probability of a control is zero) will be returned as \code{NA}.

  If \code{at = "points"} the calculation is performed
  only at the data points \eqn{x_i}{x[i]}. By default
  the result is a vector of values
  \eqn{p(x_i)}{p(x[i])} giving the estimated probability of a case
  at each data point, or a matrix of values 
  \eqn{p_j(x_i)}{p[j](x[i])} giving the estimated probability of
  each possible type \eqn{j} at each data point.
  If \code{relative=TRUE} then the relative risks
  \eqn{r(x_i)}{r(x[i])} or \eqn{r_j(x_i)}{r[j](x[i])} are
  returned.
  An infinite value of relative risk (arising because the
  probability of a control is zero) will be returned as \code{Inf}.

  Probabilities and risks are computed from the fitted intensity of the model,
  using \code{\link{predict.ppm}}.
  If \code{se=TRUE} then standard errors will also be computed,
  based on asymptotic theory, using \code{\link{vcov.ppm}}.
}
\value{
  If \code{se=FALSE} (the default), the format is described below.
  If \code{se=TRUE}, the result is a list of two entries,
  \code{estimate} and \code{SE}, each having the format described below.
  
  If \code{X} consists of only two types of points,
  and if \code{casecontrol=TRUE},
  the result is a pixel image (if \code{at="pixels"})
  or a vector (if \code{at="points"}).
  The pixel values or vector values
  are the probabilities of a case if \code{relative=FALSE},
  or the relative risk of a case (probability of a case divided by the
  probability of a control) if \code{relative=TRUE}.

  If \code{X} consists of more than two types of points,
  or if \code{casecontrol=FALSE}, the result is:
  \itemize{
    \item (if \code{at="pixels"})
    a list of pixel images, with one image for each possible type of point.
    The result also belongs to the class \code{"listof"} so that it can
    be printed and plotted.
    \item
    (if \code{at="points"})
    a matrix of probabilities, with rows corresponding to
    data points \eqn{x_i}{x[i]}, and columns corresponding
    to types \eqn{j}.
  }
  The pixel values or matrix entries
  are the probabilities of each type of point if \code{relative=FALSE},
  or the relative risk of each type (probability of each type divided by the
  probability of a control) if \code{relative=TRUE}.

  If \code{relative=FALSE}, the resulting values always lie between 0
  and 1. If \code{relative=TRUE}, the results are either non-negative
  numbers, or the values \code{Inf} or \code{NA}. 
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\seealso{
  There is another method \code{\link{relrisk.ppp}} for point pattern datasets
  which computes \emph{nonparametric} estimates of relative risk
  by kernel smoothing.

  See also
  \code{\link{relrisk}}, 
  \code{\link{relrisk.ppp}},
  \code{\link{ppm}}
}
\examples{
  fit <- ppm(chorley ~ marks * (x+y))
  rr <- relrisk(fit, relative=TRUE, control="lung", se=TRUE)
  plot(rr$estimate)
  plot(rr$SE)
}
\keyword{spatial}
\keyword{models}
