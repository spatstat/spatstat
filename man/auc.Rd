\name{auc}
\alias{auc}
\title{
  Area Under ROC Curve
}
\description{
  Compute the AUC (area under the Receiver Operating Characteristic
  curve) for a fitted point process model.
}
\usage{
auc(model)
}
\arguments{
  \item{model}{
    Fitted point process model (object of class \code{"ppm"}).
  }
}
\details{
  This command computes the AUC, the area under the Receiver Operating
  Characteristic curve, for the fitted \code{model}.
  
  The AUC is a numerical index that measures the ability of the
  fitted model to discriminate or separate the spatial domain into areas of high
  and low density of points.

  Suppose \eqn{\lambda(u)}{\lambda(u)} is the intensity function of the model,
  and \eqn{x} is the original point pattern dataset to which the
  model was fitted. Let \eqn{x_i}{x[i]} be a randomly-chosen data point,
  and \eqn{U} a randomly-selected location in the study region.
  The AUC is the probability that
  \eqn{\lambda(x_i) > \lambda(U)}{lambda(x[i]) > lambda(U)}.
  That is, AUC is the probability that a randomly-selected data point
  has higher predicted intensity than does a randomly-selected spatial
  location. The AUC is a number between 0 and 1. A value of 0.5
  indicates a complete lack of discriminatory power.
  
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
Adrian Baddeley
\email{Adrian.Baddeley@uwa.edu.au}
\url{http://www.maths.uwa.edu.au/~adrian/}
Rolf Turner
\email{r.turner@auckland.ac.nz}
and Ege Rubak
\email{rubak@math.aau.dk}
}
\examples{
  fit <- ppm(swedishpines ~ x+y)
  auc(fit)

}
\keyword{spatial}
