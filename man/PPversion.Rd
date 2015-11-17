\name{PPversion}
\alias{PPversion}
\alias{QQversion}
\title{
  Transform a Function into its P-P or Q-Q Version
}
\description{
  Given a function object \code{f} containing both the estimated
  and theoretical versions of a summary function, these operations
  combine the estimated and theoretical functions into a new function.
  When plotted, the new function gives either the P-P plot or Q-Q plot
  of the original \code{f}.
}
\usage{
PPversion(f, theo = "theo", columns = ".")

QQversion(f, theo = "theo", columns = ".")
}
\arguments{
  \item{f}{
    The function to be transformed. An object of class \code{"fv"}.
  }
  \item{theo}{
    The name of the column of \code{f} that should be treated as the
    theoretical value of the function.
  }
  \item{columns}{
    Character vector, specifying the columns of \code{f}
    to which the transformation will be applied.
    Either a vector of names of columns of \code{f},
    or one of the abbreviations recognised by \code{\link{fvnames}}.
  }
}
\details{
  The argument \code{f} should be an object of class \code{"fv"},
  containing both empirical estimates \eqn{\widehat f(r)}{fhat(r)}
  and a theoretical value \eqn{f_0(r)}{f0(r)} for a summary function.

  The \emph{P--P version} of \code{f} is the function
  \eqn{g(x) = \widehat f (f_0^{-1}(x))}{g(x) = fhat(f0^(-1)(x))}
  where \eqn{f_0^{-1}}{f0^(-1)} is the inverse function of
  \eqn{f_0}{f0}.
  A plot of \eqn{g(x)} against \eqn{x} 
  is equivalent to a plot of \eqn{\widehat f(r)}{fhat(r)} against
  \eqn{f_0(r)}{f0(r)} for all \eqn{r}.
  If \code{f} is a cumulative distribution function (such as the
  result of \code{\link{Fest}} or \code{\link{Gest}}) then
  this is a P--P plot, a plot of the observed versus theoretical
  probabilities for the distribution.
  The diagonal line \eqn{y=x}
  corresponds to perfect agreement between observed and theoretical
  distribution.

  The \emph{Q--Q version} of \code{f} is the function
  \eqn{h(x) = f_0^{-1}(\widehat f(x))}{f0^(-1)(fhat(x))}.
  If \code{f} is a cumulative distribution function,
  a plot of \eqn{h(x)} against \eqn{x}
  is a Q--Q plot, a plot of the observed versus theoretical
  quantiles of the distribution.
  The diagonal line \eqn{y=x}
  corresponds to perfect agreement between observed and theoretical
  distribution.
  Another straight line corresponds to the situation where the
  observed variable is a linear transformation of the theoretical variable.
  For a point pattern \code{X}, the Q--Q version of \code{Kest(X)} is
  essentially equivalent to \code{Lest(X)}.
}
\value{
  Another object of class \code{"fv"}.
}
\author{
  Tom Lawrence
  and Adrian Baddeley.
  
  Implemented by
  \spatstatAuthors.
}
\seealso{
  \code{\link{plot.fv}}
}
\examples{
  G <- Gest(redwoodfull)
  plot(PPversion(G))
  plot(QQversion(G))
}
\keyword{spatial}
\keyword{nonparametric}
\keyword{manip}
