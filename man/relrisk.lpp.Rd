\name{relrisk.lpp}
\alias{relrisk.lpp}
\title{
  Nonparametric Estimate of Spatially-Varying Relative Risk on a Network
}
\description{
  Given a multitype point pattern on a linear network,
  this function estimates the
  spatially-varying probability of each type of point, or the ratios of
  such probabilities, using kernel smoothing.
}
\usage{
\method{relrisk}{lpp}(X, sigma, ..., 
           at = c("pixels", "points"),
           relative=FALSE,
           adjust=1, 
           casecontrol=TRUE, control=1, case,
           finespacing=FALSE)
}
\arguments{
  \item{X}{
    A multitype point pattern (object of class \code{"lpp"}
    which has factor valued marks).
  }
  \item{sigma}{
    The numeric value of the smoothing bandwidth
    (the standard deviation of Gaussian smoothing kernel)
    passed to \code{\link{density.lpp}}.
    Alternatively \code{sigma} may be a function which can be used
    to select the bandwidth. See Details.
  }
  \item{\dots}{
    Arguments passed to \code{\link{density.lpp}} to control the
    pixel resolution.
  }
  \item{at}{
    Character string specifying whether to compute the probability values
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
  \item{adjust}{
    Optional. Adjustment factor for the bandwidth \code{sigma}.
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
  \item{finespacing}{
    Logical value specifying whether to use a finer spatial
    resolution (with longer computation time but higher accuracy).
  }
}
\details{
  The command \code{\link{relrisk}} is generic and can be used to
  estimate relative risk in different ways.
  
  This function \code{relrisk.lpp} is the method for point patterns
  on a linear network (objects of class \code{"lpp"}).
  It computes \emph{nonparametric} estimates of relative risk
  by kernel smoothing.

  If \code{X}  is a bivariate point pattern
  (a multitype point pattern consisting of two types of points)
  then by default,
  the points of the first type (the first level of \code{marks(X)})
  are treated as controls or non-events, and points of the second type
  are treated as cases or events. Then by default this command computes
  the spatially-varying \emph{probability} of a case,
  i.e. the probability \eqn{p(u)}
  that a point at location \eqn{u} on the network
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
  that a point at location \eqn{u} on the network
  will belong to type \eqn{j}.
  If \code{relative=TRUE}, the command computes the
  \emph{relative risk} of an event of type \eqn{j}
  relative to a control,
  \eqn{r_j(u) = p_j(u)/p_k(u)}{r[j](u) = p[j](u)/p[k](u)},
  where events of type \eqn{k} are treated as controls.
  The argument \code{control} determines which type \eqn{k}
  is treated as a control.

  If \code{at = "pixels"} the calculation is performed for
  every location \eqn{u} on a fine pixel grid over the network, and the result
  is a pixel image on the network representing the function \eqn{p(u)},
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

  Estimation is performed by a Nadaraja-Watson type kernel
  smoother (McSwiggan et al., 2019).

  The smoothing bandwidth \code{sigma}
  should be a single numeric value, giving the standard
  deviation of the isotropic Gaussian kernel.
  If \code{adjust} is given, the smoothing bandwidth will be
  \code{adjust * sigma} before the computation of relative risk.

  Alternatively, \code{sigma} may be a function that can be applied
  to the point pattern \code{X} to select a bandwidth; the function
  must return a single numerical value; examples include the
  functions \code{\link{bw.relrisklpp}} and \code{\link{bw.scott.iso}}.

  Accuracy depends on the spatial resolution of the density
  computations. If the arguments \code{dx} and \code{dt} are present,
  they are passed to \code{\link{density.lpp}} to determine the
  spatial resolution. Otherwise, the spatial resolution is determined
  by a default rule that depends on \code{finespacing} and \code{sigma}.
  If \code{finespacing=FALSE} (the default), the spatial resolution is
  equal to the default resolution for pixel images.
  If \code{finespacing=TRUE}, the spatial resolution is much finer
  and is determined by a rule which guarantees higher accuracy,
  but takes a longer time.
}
\value{
   If \code{X} consists of only two types of points,
  and if \code{casecontrol=TRUE},
  the result is a pixel image on the network (if \code{at="pixels"})
  or a vector (if \code{at="points"}).
  The pixel values or vector values
  are the probabilities of a case if \code{relative=FALSE},
  or the relative risk of a case (probability of a case divided by the
  probability of a control) if \code{relative=TRUE}.

  If \code{X} consists of more than two types of points,
  or if \code{casecontrol=FALSE}, the result is:
  \itemize{
    \item (if \code{at="pixels"})
    a list of pixel images on the network,
    with one image for each possible type of point.
    The result also belongs to the class \code{"solist"} so that it can
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
\seealso{
  \code{\link{relrisk}}
}
\examples{
   ## case-control data: 2 types of points
   set.seed(2020)
   X <- superimpose(A=runiflpp(20, simplenet),
                    B=runifpointOnLines(20, as.psp(simplenet)[1]))
   plot(X)
   plot(relrisk(X, 0.2))
   plot(relrisk(X, 0.2, case="B"))
   head(relrisk(X, 0.2, at="points"))
   ## cross-validated bandwidth selection
   plot(relrisk(X, bw.relrisklpp, hmax=0.3))

   ## more than 2 types
   if(interactive()) {
     U <- chicago
     sig <- 170
   } else {
     U <- do.call(superimpose,
                  split(chicago)[c("theft", "cartheft", "burglary")])
     sig <- 50
   }
   plot(relrisk(U, sig))
   head(relrisk(U, sig, at="points"))
   plot(relrisk(U, sig, relative=TRUE, control="theft"))
}
\references{
  McSwiggan, G., Baddeley, A. and Nair, G. (2019)
  Estimation of relative risk for events on a linear network.
  \emph{Statistics and Computing} \bold{30} (2) 469--484.
}
\author{
  Greg McSwiggan and \adrian.
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
