\name{plot.slrm}
\Rdversion{1.1}
\alias{plot.slrm}
\title{
  Plot a Fitted Spatial Logistic Regression
}
\description{
  Plots a fitted Spatial Logistic Regression model.
}
\usage{
  \method{plot}{slrm}(x, ..., type = "intensity")
}

\arguments{
  \item{x}{
    a fitted spatial logistic regression model.
    An object of class \code{"slrm"}.
  }
  \item{\dots}{
    Extra arguments passed to \code{\link{plot.im}} to
    control the appearance of the plot.
  }
  \item{type}{
    Character string (partially) matching one of 
    \code{"probabilities"}, \code{"intensity"} or \code{"link"}.
  }
}
\details{
  This is a method for \code{\link{plot}} for fitted spatial logistic
  regression models (objects of class \code{"slrm"}, usually obtained
  from the function \code{\link{slrm}}).

  This function plots the result of \code{\link{predict.slrm}}.
}
\value{
  None.
}
\seealso{
  \code{\link{slrm}}, 
  \code{\link{predict.slrm}}, 
  \code{\link{plot.im}}
}
\examples{
   data(copper)
   X <- copper$SouthPoints
   Y <- copper$SouthLines
   Z <- distmap(Y)
   fit <- slrm(X ~ Z)
   plot(fit)
   plot(fit, type="link")
}
\author{\adrian
  \email{adrian@maths.uwa.edu.au}
  
  and \rolf
  
}
\keyword{spatial}
\keyword{hplot}
\keyword{models}
