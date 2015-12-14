\name{markcrosscorr}
\alias{markcrosscorr}
\title{
  Mark Cross-Correlation Function
}
\description{
  Given a spatial point pattern with several columns of marks,
  this function computes the mark correlation function between
  each pair of columns of marks.
}
\usage{
  markcrosscorr(X, r = NULL,
                correction = c("isotropic", "Ripley", "translate"),
                method = "density", \dots, normalise = TRUE, Xname = NULL)
}
\arguments{
  \item{X}{The observed point pattern.
    An object of class \code{"ppp"} or something acceptable to
    \code{\link{as.ppp}}. 
  }
  \item{r}{Optional. Numeric vector. The values of the argument \eqn{r}
    at which the mark correlation function 
    \eqn{k_f(r)}{k[f](r)} should be evaluated.
    There is a sensible default.
  }
  \item{correction}{
    A character vector containing any selection of the
    options \code{"isotropic"}, \code{"Ripley"}, \code{"translate"},
    \code{"translation"}, \code{"none"} or \code{"best"}.
    It specifies the edge correction(s) to be applied.
  }
  \item{method}{
    A character vector indicating the user's choice of
    density estimation technique to be used. Options are
    \code{"density"}, 
    \code{"loess"},
    \code{"sm"} and \code{"smrep"}.
  }
  \item{\dots}{
    Arguments passed to the density estimation routine
    (\code{\link{density}}, \code{\link{loess}} or \code{sm.density})
    selected by \code{method}.
  }
  \item{normalise}{
    If \code{normalise=FALSE},
    compute only the numerator of the expression for the
    mark correlation.
  }
  \item{Xname}{
    Optional character string name for the dataset \code{X}.
  }
}
\details{
  First, all columns of marks are converted to numerical values.
  A factor with \eqn{m} possible levels is converted to
  \eqn{m} columns of dummy (indicator) values.

  Next, each pair of columns is considered, and the mark
  cross-correlation is defined as
  \deqn{
    k_{mm}(r) = \frac{E_{0u}[M_i(0) M_j(u)]}{E[M_i,M_j]}
  }{
    k[mm](r) = E[0u](M(i,0) * M(j,u))/E(Mi * Mj)
  }
  where \eqn{E_{0u}}{E[0u]} denotes the conditional expectation
  given that there are points of the process at the locations
  \eqn{0} and \eqn{u} separated by a distance \eqn{r}.
  On the numerator,
  \eqn{M_i(0)}{M(i,0)} and \eqn{M_j(u)}{M(j,u)}
  are the marks attached to locations \eqn{0} and \eqn{u} respectively
  in the \eqn{i}th and \eqn{j}th columns of marks respectively.
  On the denominator, \eqn{M_i}{Mi} and \eqn{M_j}{Mj} are
  independent random values drawn from the
  \eqn{i}th and \eqn{j}th columns of marks, respectively,
  and \eqn{E} is the usual expectation.
  
  Note that \eqn{k_{mm}(r)}{k[mm](r)} is not a ``correlation''
  in the usual statistical sense. It can take any 
  nonnegative real value. The value 1 suggests ``lack of correlation'':
  if the marks attached to the points of \code{X} are independent
  and identically distributed, then
  \eqn{k_{mm}(r) \equiv 1}{k[mm](r) =  1}.

  The argument \code{X} must be a point pattern (object of class
  \code{"ppp"}) or any data that are acceptable to \code{\link{as.ppp}}.
  It must be a marked point pattern.

  The cross-correlations are estimated in the same manner as
  for \code{\link{markcorr}}.
}
\value{
  A function array (object of class \code{"fasp"}) containing
  the mark cross-correlation functions for each possible pair
  of columns of marks.
}
\author{
  Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  Rolf Turner \email{r.turner@auckland.ac.nz}
  
  and Ege Rubak \email{rubak@math.aau.dk}
  
}
\seealso{
  \code{\link{markcorr}}
}
\examples{
  # The dataset 'betacells' has two columns of marks:
  #       'type' (factor)
  #       'area' (numeric)
  if(interactive()) plot(betacells)
  plot(markcrosscorr(betacells))
}
\keyword{spatial}
\keyword{nonparametric}
