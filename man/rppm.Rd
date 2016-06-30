\name{rppm}
\alias{rppm}
\title{
  Recursively Partitioned Point Process Model
}
\description{
  Fits a recursive partition model to point pattern data.
}
\usage{
  rppm(\dots, rpargs=list())
}
\arguments{
  \item{\dots}{
    Arguments passed to \code{\link{ppm}}
    specifying the point pattern data and the explanatory covariates.
  }
  \item{rpargs}{
    Optional list of arguments passed to \code{\link[rpart]{rpart}}
    controlling the recursive partitioning procedure.
  }
}
\details{
  This function attempts to find a simple rule for predicting
  low and high intensity regions of points in a point pattern,
  using explanatory covariates.
  
  The arguments \code{\dots} specify the point pattern data
  and explanatory covariates in the same way as they would be
  in the function \code{\link{ppm}}. 

  The recursive partitioning algorithm \code{\link[rpart]{rpart}}
  is then used to find a partitioning rule.
}
\value{
  An object of class \code{"rppm"}. There are methods
  for \code{print}, \code{plot}, \code{fitted}, \code{predict} and
  \code{prune} for this class.
}
\references{
  Breiman, L., Friedman, J. H., Olshen, R. A., and Stone, C. J. (1984)
  \emph{Classification and Regression Trees}. Wadsworth.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{plot.rppm}}, 
  \code{\link{predict.rppm}}, 
  \code{\link{prune.rppm}}.
}
\examples{
    # New Zealand trees data: trees planted along border
    # Use covariates 'x', 'y'
    nzfit <- rppm(nztrees ~ x + y)
    nzfit
    prune(nzfit, cp=0.035)
    # Murchison gold data: numeric and logical covariates
    mur <- solapply(murchison, rescale, s=1000, unitname="km")
    mur$dfault <- distfun(mur$faults)
    # 
    mfit <- rppm(gold ~ dfault + greenstone, data=mur)
    mfit
    # Gorillas data: factor covariates
    #          (symbol '.' indicates 'all variables')
    gfit <- rppm(unmark(gorillas) ~ . , data=gorillas.extra)
    gfit
}
\keyword{spatial}
\keyword{models}
