\name{relrisk}
\alias{relrisk}
\title{
  Nonparametric Estimate of Spatially-Varying Relative Risk
}
\description{
  Given a multitype point pattern, this function estimates the
  spatially-varying probability of each type of point, using
  kernel smoothing. The default smoothing bandwidth is selected by
  cross-validation.
}
\usage{
relrisk(X, sigma = NULL, ..., varcov = NULL, at = "pixels", casecontrol=TRUE)
}
\arguments{
  \item{X}{
    A multitype point pattern (object of class \code{"ppp"}
    which has factor valued marks).
  }
  \item{sigma}{
    Optional. Standard deviation of isotropic Gaussian smoothing kernel.
  }
  \item{\dots}{
    Arguments passed to \code{\link{bw.relrisk}} to select the
    bandwidth, or passed to \code{\link{density.ppp}} to control the
    pixel resolution.
  }
  \item{varcov}{
    Optional. Variance-covariance matrix of anisotopic Gaussian
    smoothing kernel. Incompatible with \code{sigma}.
  }
  \item{at}{
    String specifying whether to compute the probability values
    at a grid of pixel locations (\code{at="pixels"}) or
    only at the points of \code{X} (\code{at="points"}).
  }
  \item{casecontrol}{
    Logical. Whether to treat a bivariate point pattern
    as consisting of cases and controls. See Details.
  }
}
\details{
  If \code{X}  is a bivariate point pattern
  (a multitype point pattern consisting of two types of points)
  then by default,
  the points of the first type (the first level of \code{marks(X)})
  are treated as controls or non-events, and points of the second type
  are treated as cases or events. Then this command computes
  the spatially-varying risk of an event,
  i.e. the probability \eqn{p(u)}
  that a point at spatial location \eqn{u}
  will be a case.

  If \code{X} is a multitype point pattern with \eqn{m > 2} types,
  or if \code{X} is a bivariate point pattern
  and \code{casecontrol=FALSE},
  then this command computes, for each type \eqn{j},
  a nonparametric estimate of
  the spatially-varying risk of an event of type \eqn{j}.
  This is the probability \eqn{p_j(u)}{p[j](u)}
  that a point at spatial location \eqn{u}
  will belong to type \eqn{j}.

  If \code{at = "pixels"} the calculation is performed for
  every spatial location \eqn{u} on a fine pixel grid, and the result
  is a pixel image representing the function \eqn{p(u)}
  or a list of pixel images representing the functions 
  \eqn{p_j(u)}{p[j](u)} for \eqn{j = 1,\ldots,m}{j = 1,...,m}.

  If \code{at = "points"} the calculation is performed
  only at the data points \eqn{x_i}{x[i]}. The result is a vector of values
  \eqn{p(x_i)}{p(x[i])} giving the estimated probability of a case
  at each data point, or a matrix of values 
  \eqn{p_j(x_i)}{p[j](x[i])} giving the estimated probability of
  each possible type \eqn{j} at each data point.

  Estimation is performed by a simple Nadaraja-Watson type kernel
  smoother (Diggle, 2003). If \code{sigma} and \code{varcov}
  are both missing or null, then the smoothing bandwidth \code{sigma}
  is selected by cross-validation using \code{\link{bw.relrisk}}.
}
\value{
  If \code{X} consists of only two types of points,
  the result is a pixel image (if \code{at="pixels"})
  or a vector of probabilities (if \code{at="points"}).

  If \code{X} consists of more than two types of points,
  the result is:
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
}
\seealso{
 \code{\link{bw.relrisk}},
 \code{\link{density.ppp}},
 \code{\link{Smooth.ppp}},
 \code{\link{eval.im}}.
 \code{\link{which.max.im}}.
}
\examples{
   data(urkiola)
   p <- relrisk(urkiola, 20)
   if(interactive()) {
      plot(p, main="proportion of oak")
      plot(eval.im(p > 0.3), main="More than 30 percent oak")
      data(lansing)
      z <- relrisk(lansing)
      plot(z, main="Lansing Woods")
      plot(which.max.im(z), main="Most common species")
   }
}
\references{
  Diggle, P.J. (2003)
  \emph{Statistical analysis of spatial point patterns},
  Second edition. Arnold.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
