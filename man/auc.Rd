\name{auc}
\alias{auc}
\alias{auc.ppp}
\alias{auc.ppm}
\title{
  Area Under ROC Curve
}
\description{
  Compute the AUC (area under the Receiver Operating Characteristic
  curve) for a fitted point process model.
}
\usage{
auc(X, \dots)

\method{auc}{ppp}(X, covariate, \dots, high = TRUE)

\method{auc}{ppm}(X, \dots)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"})
    or fitted point process model (object of class \code{"ppm"}).
  }
  \item{covariate}{
    Spatial covariate. Either a \code{function(x,y)},
    a pixel image (object of class \code{"im"}), or
    one of the strings \code{"x"} or \code{"y"} indicating the
    Cartesian coordinates.    
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} controlling the
    pixel resolution for calculations.
  }
  \item{high}{
    Logical value indicating whether the threshold operation
    should favour high or low values of the covariate.
  }
}
\details{
  This command computes the AUC, the area under the Receiver Operating
  Characteristic curve. The ROC itself is computed by \code{\link{roc}}.

  For a point pattern \code{X} and a covariate \code{Z}, the
  AUC is a numerical index that measures the ability of the 
  covariate to separate the spatial domain
  into areas of high and low density of points.
  Let \eqn{x_i}{x[i]} be a randomly-chosen data point from \code{X}
  and \eqn{U} a randomly-selected location in the study region.
  The AUC is the probability that
  \eqn{Z(x_i) > Z(U)}{Z(x[i]) > Z(U)}
  assuming \code{high=TRUE}.
  That is, AUC is the probability that a randomly-selected data point
  has a higher value of the covariate \code{Z} than does a
  randomly-selected spatial location. The AUC is a number between 0 and 1. 
  A value of 0.5 indicates a complete lack of discriminatory power.
  
  For a fitted point process model \code{X},
  the AUC measures the ability of the
  fitted model intensity to separate the spatial domain
  into areas of high and low density of points.
  Suppose \eqn{\lambda(u)}{\lambda(u)} is the intensity function of the model.
  The AUC is the probability that
  \eqn{\lambda(x_i) > \lambda(U)}{lambda(x[i]) > lambda(U)}.
  That is, AUC is the probability that a randomly-selected data point
  has higher predicted intensity than does a randomly-selected spatial
  location.
  The AUC is \bold{not} a measure of the goodness-of-fit of the model
  (Lobo et al, 2007).
}
\value{
  A numeric vector of length 2 giving the AUC value 
  and the theoretically expected AUC value for this model.
}
\references{
  Lobo, J.M.,
  \ifelse{latex}{\out{Jim{\'e}nez}}{Jimenez}-Valverde, A.
  and Real, R. (2007)
  AUC: a misleading measure of the performance of predictive
  distribution models.
  \emph{Global Ecology and Biogeography} \bold{17}(2) 145--151.

  Nam, B.-H. and D'Agostino, R. (2002)
  Discrimination index, the area under the {ROC} curve.
  Pages 267--279 in 
  Huber-Carol, C., Balakrishnan, N., Nikulin, M.S. 
  and Mesbah, M., \emph{Goodness-of-fit tests and model validity},
  \ifelse{latex}{\out{Birkh{\"a}user}}{Birkhauser}, Basel.
}
\author{
\adrian


\rolf

and \ege

}
\seealso{
  \code{\link{roc}}
}
\examples{
  fit <- ppm(swedishpines ~ x+y)
  auc(fit)
  auc(swedishpines, "x")
}
\keyword{spatial}
