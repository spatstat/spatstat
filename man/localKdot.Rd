\name{localKdot}
\alias{localKdot}
\alias{localLdot}
\title{Local Multitype K Function (Dot-Type)}
\description{
  for a multitype point pattern,
  computes the dot-type version of
  the local K function.
}
\usage{
  localKdot(X, from, \dots, rmax = NULL,
              correction = "Ripley", verbose = TRUE, rvalue=NULL)
  localLdot(X, from, \dots, rmax = NULL, correction = "Ripley")
}
\arguments{
  \item{X}{
    A multitype point pattern (object of class \code{"ppp"}
    with marks which are a factor).
  }
  \item{\dots}{
    Further arguments passed from \code{localLdot} to
    \code{localKdot}.
  }
  \item{rmax}{
    Optional. Maximum desired value of the argument \eqn{r}.
  }
  \item{from}{
    Type of points from which distances should be measured.
    A single value;
    one of the possible levels of \code{marks(X)},
    or an integer indicating which level.
  }
  \item{correction}{
    String specifying the edge correction to be applied.
    Options are \code{"none"}, \code{"translate"}, \code{"translation"},
    \code{"Ripley"},
    \code{"isotropic"} or \code{"best"}.
    Only one correction may be specified.
  }
  \item{verbose}{
    Logical flag indicating whether to print progress
    reports during the calculation.
  }
  \item{rvalue}{Optional. A \emph{single} value of the distance argument
    \eqn{r} at which the function L or K should be computed.
  }
}
\details{
  Given a multitype spatial point pattern \code{X},
  the local dot-type \eqn{K} function \code{localKdot}
  is the local version of the multitype \eqn{K} function
  \code{\link{Kdot}}. 
  Recall that \code{Kdot(X, from)} is a sum of contributions
  from all pairs of points in \code{X} where
  the first point belongs to \code{from}.
  The \emph{local} dot-type \eqn{K}
  function is defined for each point \code{X[i]} that belongs to
  type \code{from}, and it consists of all the contributions to
  the dot-type \eqn{K} function that originate from point \code{X[i]}:
  \deqn{
    K_{i,from,to}(r) = \sqrt{\frac a {(n-1) \pi} \sum_j e_{ij}}
  }{
    K[i,from,to](r) = sqrt( (a/((n-1)* pi)) * sum[j] e[i,j])
  }
  where the sum is over all points \eqn{j \neq i}{j != i}
  that lie within a distance \eqn{r} of the \eqn{i}th point, 
  \eqn{a} is the area of the observation window, \eqn{n} is the number
  of points in \code{X}, and \eqn{e_{ij}}{e[i,j]} is an edge correction
  term (as described in \code{\link{Kest}}).
  The value of \eqn{K_{i,from}(r)}{K[i,from](r)}
  can also be interpreted as one
  of the summands that contributes to the global estimate of the
  \code{\link{Kdot}} function.

  By default, the function \eqn{K_{i,from}(r)}{K[i,from](r)}
  is computed for a range of \eqn{r} values
  for each point \eqn{i} belonging to type \code{from}.
  The results are stored as a function value
  table (object of class \code{"fv"}) with a column of the table
  containing the function estimates for each point of the pattern
  \code{X} belonging to type \code{from}.

  Alternatively, if the argument \code{rvalue} is given, and it is a
  single number, then the function will only be computed for this value
  of \eqn{r}, and the results will be returned as a numeric vector,
  with one entry of the vector for each point of the pattern \code{X}
  belonging to type \code{from}.

  The local dot-type \eqn{L} function \code{localLdot}
  is computed by applying the transformation
  \eqn{L(r) = \sqrt{K(r)/(2\pi)}}{L(r) = sqrt(K(r)/(2*pi))}.
}
\value{
  If \code{rvalue} is given, the result is a numeric vector
  of length equal to the number of points in the point pattern
  that belong to type \code{from}.

  If \code{rvalue} is absent, the result is 
  an object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.
  Essentially a data frame containing columns
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the function \eqn{K} has been  estimated
  }
  \item{theo}{the theoretical value \eqn{K(r) = \pi r^2}{K(r) = pi * r^2}
    or \eqn{L(r)=r} for a stationary Poisson process
  }
  together with columns containing the values of the
  neighbourhood density function for each point in the pattern.
  Column \code{i} corresponds to the \code{i}th point
  of type \code{from}.
  The last two columns contain the \code{r} and \code{theo} values.
}
\seealso{
  \code{\link{Kdot}},
  \code{\link{Ldot}},
  \code{\link{localK}},
  \code{\link{localL}}.
}
\examples{
  X <- amacrine

  # compute all the local Ldot functions
  L <- localLdot(X)

  # plot all the local Ldot functions against r
  plot(L, main="local Ldot functions for amacrine", legend=FALSE)

  # plot only the local L function for point number 7
  plot(L, iso007 ~ r)
  
  # compute the values of L(r) for r = 0.1 metres
  L12 <- localLdot(X, rvalue=0.1)
}
\author{
  \ege and \adrian.
}
\keyword{spatial}
\keyword{nonparametric}

