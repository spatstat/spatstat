\name{ppmInfluence}
\alias{ppmInfluence}
\title{
  Leverage and Influence Measures for Spatial Point Process Model
}
\description{
  Calculates all the leverage and
  influence measures described in \code{\link{influence.ppm}},
  \code{\link{leverage.ppm}} and \code{\link{dfbetas.ppm}}.
}
\usage{
   ppmInfluence(fit,
                what = c("leverage", "influence", "dfbetas"),
                \dots,
                iScore = NULL, iHessian = NULL, iArgs = NULL,
                drop = FALSE,
                fitname = NULL)
}
\arguments{
  \item{fit}{
    A fitted point process model of class \code{"ppm"}.
  }
  \item{what}{
    Character vector specifying which quantities are to be calculated.
    Default is to calculate all quantities.
  }
  \item{\dots}{
    Ignored.
  }
  \item{iScore,iHessian}{
    Components of the score vector and Hessian matrix for
    the irregular parameters, if required. See Details.
  }
  \item{iArgs}{
    List of extra arguments for the functions \code{iScore},
    \code{iHessian} if required.
  }
  \item{drop}{
    Logical. Whether to include (\code{drop=FALSE}) or
    exclude (\code{drop=TRUE}) contributions from quadrature
    points that were not used to fit the model.
  }
  \item{fitname}{
    Optional character string name for the fitted model \code{fit}.
  }
}
\details{
  This function calculates all the
  leverage and influence measures
  described in \code{\link{influence.ppm}}, \code{\link{leverage.ppm}}
  and \code{\link{dfbetas.ppm}}.

  When analysing large datasets, the user can
  call \code{ppmInfluence} to perform the calculations efficiently,
  then extract the leverage and influence values as desired.
  For example the leverage can be extracted either as
  \code{result$leverage} or \code{leverage(result)}.
  
  If the point process model trend has irregular parameters that were
  fitted (using \code{\link{ippm}})
  then the influence calculation requires the first and second
  derivatives of the log trend with respect to the irregular parameters. 
  The argument \code{iScore} should be a list,
  with one entry for each irregular parameter, of \R functions that compute the
  partial derivatives of the log trend (i.e. log intensity or
  log conditional intensity) with respect to each irregular
  parameter. The argument \code{iHessian} should be a list,
  with \eqn{p^2} entries where \eqn{p} is the number of irregular
  parameters, of \R functions that compute the second order
  partial derivatives of the
  log trend with respect to each pair of irregular parameters.
}
\value{
  A list containing the leverage and influence measures specified by
  \code{what}. The result also belongs to the class \code{"ppmInfluence"}. 
}
\author{
  \adrian
}
\seealso{
  \code{\link{leverage.ppm}},
  \code{\link{influence.ppm}},
  \code{\link{dfbetas.ppm}}
}
\examples{
   X <- rpoispp(function(x,y) { exp(3+3*x) })
   fit <- ppm(X ~ x+y)
   fI <- ppmInfluence(fit)
   fI$influence
   influence(fI) 
   fI$leverage
   fI$dfbetas
}
\keyword{spatial}
\keyword{models}
