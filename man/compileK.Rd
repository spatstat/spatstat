\name{compileK}
\alias{compileK}
\alias{compilepcf}
\title{
  Generic Calculation of K Function and Pair Correlation Function
}
\description{
  Low-level functions which 
  calculate the estimated \eqn{K} function
  and estimated pair correlation function
  (or any similar functions)
  from a matrix of pairwise distances and optional weights.
}
\usage{
compileK(D, r, weights = NULL, denom = 1,
         check = TRUE, ratio = FALSE, fname = "K")

compilepcf(D, r, weights = NULL, denom = 1,
         check = TRUE, endcorrect = TRUE, \dots, fname = "g")
}
\arguments{
  \item{D}{
    A square matrix giving the distances between all pairs of points.
  }
  \item{r}{
    An equally spaced, finely spaced sequence of distance values.
  }
  \item{weights}{
    Optional numerical weights for the pairwise distances.
    A numeric matrix with the same dimensions as \code{D}.
    If absent, the weights are taken to equal 1.
  }
  \item{denom}{
    Denominator for the estimator.
    A single number, or a numeric vector with the same length
    as \code{r}. See Details.
  }
  \item{check}{
    Logical value specifying whether to check that \code{D} is a
    valid matrix of pairwise distances.
  }
  \item{ratio}{
    Logical value indicating whether to store ratio information.
    See Details.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link[stats]{density.default}}
    controlling the kernel smoothing.
  }
  \item{endcorrect}{
    Logical value indicating whether to apply End Correction of
    the pair correlation estimate at \code{r=0}.
  }
  \item{fname}{
    Character string giving the name of the function being estimated.
  }
}
\details{
  These low-level functions construct estimates of
  the \eqn{K} function or pair correlation function,
  or any similar functions, given only the matrix of pairwise
  distances and optional weights associated with these distances.

  These functions are useful for code development and for teaching,
  because they perform a common task, and do the housekeeping required to
  make an object of class \code{"fv"} that represents the estimated
  function. However, they are not very efficient. 
  
  \code{compileK} calculates the weighted estimate
  of the \eqn{K} function,
  \deqn{
    \hat K(r) = (1/v(r)) \sum_i \sum_j 1\{ d_{ij} \le r\} w_{ij}
  }{
    K(r) = (1/v(r)) sum[i] sum[j] 1(d[i,j] <= r) w[i,j]
  }
  and \code{compilepcf} calculates the weighted estimate of the
  pair correlation function,
  \deqn{
    \hat g(r) = (1/v(r)) \sum_i \sum_j \kappa( d_{ij} - r ) w_{ij}
  }{
    g(r) = (1/v(r)) sum[i] sum[j] kappa ( d[i,j] - r) w[i,j]
  }
  where \eqn{d_{ij}}{d[i,j]} is the distance between spatial points
  \eqn{i} and \eqn{j}, with corresponding weight \eqn{w_{ij}}{w[i,j]},
  and \eqn{v(r)} is a specified denominator. Here \eqn{\kappa}{kappa}
  is a fixed-bandwidth smoothing kernel. 

  For a point pattern in two dimensions, the usual denominator \eqn{v(r)}
  is constant for the \eqn{K} function, and proportional to \eqn{r}
  for the pair correlation function. See the Examples.

  The result is an object of class \code{"fv"} representing the
  estimated function. This object has only one column of function
  values. Additional columns (such as a column giving the theoretical
  value) must be added by the user, with the aid of
  \code{\link{bind.fv}}.
}
\value{
  An object of class \code{"fv"} representing the estimated function.
}
\author{
  \adrian
}
\seealso{
  \code{\link{Kest}},
  \code{\link{pcf}} for definitions of the \eqn{K} function
  and pair correlation function.
  
  \code{\link{bind.fv}} to add more columns.
}
\examples{
  X <- japanesepines
  D <- pairdist(X)
  Wt <- edge.Ripley(X, D)
  lambda <- intensity(X)
  a <- (npoints(X)-1) * lambda
  r <- seq(0, 0.25, by=0.01)
  K <- compileK(D=D, r=r, weights=Wt, denom=a)
  g <- compilepcf(D=D, r=r, weights=Wt, denom= a * 2 * pi * r)
}
\keyword{spatial}
\keyword{nonparametric}
