\name{anova.slrm}
\Rdversion{1.1}
\alias{anova.slrm}
\title{
  Analysis of Deviance for Spatial Logistic Regression Models
}
\description{
  Performs Analysis of Deviance for two or more fitted Spatial Logistic
  Regression models.
}
\usage{
  \method{anova}{slrm}(object, ..., test = NULL)
}
\arguments{
  \item{object}{
    a fitted spatial logistic regression model.
    An object of class \code{"slrm"}.
  }
  \item{\dots}{
    additional objects of the same type (optional).
  }
  \item{test}{
    a character string, (partially) matching one of
    \code{"Chisq"}, \code{"F"} or \code{"Cp"}, indicating the
    reference distribution that should be used to compute
    \eqn{p}-values.
  }
}
\details{
  This is a method for \code{\link[stats]{anova}} for fitted spatial logistic
  regression models (objects of class \code{"slrm"}, usually obtained
  from the function \code{\link{slrm}}).

  The output shows the deviance differences (i.e. 2 times log
  likelihood ratio), the difference in degrees of freedom, and (if
  \code{test="Chi"}) the two-sided \eqn{p}-values for the chi-squared tests.
  Their interpretation is very similar to that
  in \code{\link[stats]{anova.glm}}.
}
\value{
  An object of class \code{"anova"}, inheriting from
  class \code{"data.frame"}, representing the analysis of deviance table.
}
\seealso{
  \code{\link{slrm}}
}
\examples{
  X <- rpoispp(42)
  fit0 <- slrm(X ~ 1)
  fit1 <- slrm(X ~ x+y)
  anova(fit0, fit1, test="Chi")
}
\author{\adrian
  \email{adrian@maths.uwa.edu.au}
  
  and \rolf
  
}
\keyword{spatial}
\keyword{models}
\keyword{methods}
