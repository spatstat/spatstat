\name{roc}
\alias{roc}
\alias{roc.ppp}
\alias{roc.ppm}
\title{
  Receiver Operating Characteristic
}
\description{
  Computes the Receiver Operating Characteristic curve
  for a point pattern or a fitted point process model.
}
\usage{
 roc(X, \dots)

\method{roc}{ppp}(X, covariate, \dots, high = TRUE)

\method{roc}{ppm}(X, \dots)
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
  This command computes Receiver Operating
  Characteristic curve. The area under the ROC is computed by \code{\link{auc}}.

  For a point pattern \code{X} and a covariate \code{Z}, the
  ROC is a plot showing the ability of the 
  covariate to separate the spatial domain
  into areas of high and low density of points.
  For each possible threshold \eqn{z}, the algorithm calculates
  the fraction \eqn{a(z)} of area in the study region where the
  covariate takes a value greater than \eqn{z}, and the
  fraction \eqn{b(z)} of data points for which the covariate value
  is greater than \eqn{z}. The ROC is a plot of \eqn{b(z)} against
  \eqn{a(z)} for all thresholds \eqn{z}. 
  
  For a fitted point process model, 
  the ROC shows the ability of the
  fitted model intensity to separate the spatial domain
  into areas of high and low density of points.
  The ROC is \bold{not} a diagnostic for the goodness-of-fit of the model
  (Lobo et al, 2007).
}
\value{
  Function value table (object of class \code{"fv"})
  which can be plotted to show the ROC curve.
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
Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}


Rolf Turner \email{r.turner@auckland.ac.nz}

and Ege Rubak \email{rubak@math.aau.dk}

}
\seealso{
\code{\link{auc}}
}
\examples{
  plot(roc(swedishpines, "x"))
  fit <- ppm(swedishpines ~ x+y)
  plot(roc(fit))
}
\keyword{spatial}
